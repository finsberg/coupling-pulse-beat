import dolfin
import pulse
import logging
import beat
import numpy as np
from dataclasses import dataclass
from typing import Callable
import numpy.typing as npt

logger = logging.getLogger(__name__)


def enlist(obj):
    try:
        return list(obj)
    except TypeError:
        return [obj]


@dataclass
class ODESolver:  # (beat.odesolver.DolfinODESolver):
    Ta: dolfin.Function
    init_states: npt.NDArray
    parameters: npt.NDArray
    fun: Callable
    num_states: int
    fun_monitor: Callable
    Ta_index: int = 0
    missing_variables: npt.NDArray | None = None
    num_missing_variables: int = 0

    def __post_init__(self):
        if np.shape(self.init_states) == self.shape:
            self._values = np.copy(self.init_states)
        else:
            self._values = np.zeros(self.shape)
            self._values.T[:] = self.init_states

        # if self.missing_variables is not None:
        #     self._missing_variables = self.missing_variables
        #     # if np.shape(self.missing_variables) == self.shape_missing_values:
        #     #     self._missing_variables = np.copy(self.missing_variables)
        #     # else:
        #     #     self._missing_variables = np.zeros(
        #     #         (len(self.missing_variables), self.num_points)
        #     #     )
        #     #     self._missing_variables.T[:] = self.missing_variables
        # else:
        #     self._missing_variables = None

        self._ode = beat.odesolver.ODESystemSolver(
            fun=self.fun,
            states=self._values,
            parameters=self.parameters,
            missing_variables=self.missing_variables,
        )
        # self._initialize_metadata()

    @property
    def shape_missing_values(self) -> tuple[int, int]:
        return (self.num_missing_variables, self.num_points)

    @property
    def shape(self) -> tuple[int, int]:
        return (self.num_states, self.num_points)

    @property
    def num_points(self) -> int:
        return self.Ta.vector().local_size()

    def step(self, t0: float, dt: float) -> None:
        self._ode.step(t0, t0 + dt)
        self.update_Ta(t0 + dt)

    def update_Ta(self, t):
        monitor = self.fun_monitor(
            t=t,
            states=self._values,
            parameters=self.parameters,
            missing_variables=self.missing_variables,
        )
        print("Ta: ", monitor[self.Ta_index])
        print(f"t = {t}")
        print("Missing variables: ", self.missing_variables)
        self.Ta.vector().set_local(monitor[self.Ta_index])


class NonlinearProblem(dolfin.NonlinearProblem):
    def __init__(
        self,
        J,
        F,
        bcs,
        output_matrix=False,
        output_matrix_path="output",
        **kwargs,
    ):
        super().__init__(**kwargs)
        self._J = J
        self._F = F

        self.bcs = enlist(bcs)
        self.output_matrix = output_matrix
        self.output_matrix_path = output_matrix_path
        self.verbose = True
        self.n = 0

    def F(self, b: dolfin.PETScVector, x: dolfin.PETScVector):
        dolfin.assemble(self._F, tensor=b)
        for bc in self.bcs:
            bc.apply(b, x)

    def J(self, A: dolfin.PETScMatrix, x: dolfin.PETScVector):
        dolfin.assemble(self._J, tensor=A)
        for bc in self.bcs:
            bc.apply(A)


class NewtonSolver(dolfin.NewtonSolver):
    def __init__(
        self,
        problem: pulse.NonlinearProblem,
        state: dolfin.Function,
        ode: ODESolver,
        update_cb=None,
        parameters=None,
    ):
        print(f"Initialize NewtonSolver with parameters: {parameters!r}")
        dolfin.PETScOptions.clear()
        self.ode = ode
        self._problem = problem
        self._state = state
        self._update_cb = update_cb
        self._prev_state = dolfin.Vector(state.vector().copy())
        self._diff = dolfin.Vector(state.vector().copy())

        # Initializing Newton solver (parent class)
        self.petsc_solver = dolfin.PETScKrylovSolver()
        super().__init__(
            self._state.function_space().mesh().mpi_comm(),
            self.petsc_solver,
            dolfin.PETScFactory.instance(),
        )
        self._handle_parameters(parameters)

    def _handle_parameters(self, parameters):
        # Setting default parameters
        params = type(self).default_solver_parameters()
        params.update(parameters)

        for k, v in params.items():
            if self.parameters.has_parameter(k):
                self.parameters[k] = v
            if self.parameters.has_parameter_set(k):
                for subk, subv in params[k].items():
                    self.parameters[k][subk] = subv
        petsc = params.pop("petsc")
        for k, v in petsc.items():
            if v is not None:
                dolfin.PETScOptions.set(k, v)
        self.newton_verbose = params.pop("newton_verbose", False)
        self.ksp_verbose = params.pop("ksp_verbose", False)
        self.debug = params.pop("debug", False)
        if self.newton_verbose:
            dolfin.set_log_level(dolfin.LogLevel.INFO)
            self.parameters["report"] = True
        if self.ksp_verbose:
            self.parameters["lu_solver"]["report"] = True
            self.parameters["lu_solver"]["verbose"] = True
            self.parameters["krylov_solver"]["monitor_convergence"] = True
            dolfin.PETScOptions.set("ksp_monitor_true_residual")
        self.linear_solver().set_from_options()
        self._residual_index = 0
        self._residuals = []
        self.parameters["convergence_criterion"] = "incremental"

    # def register_datacollector(self, datacollector):
    #     self._datacollector = datacollector

    @staticmethod
    def default_solver_parameters():
        return {
            "petsc": {
                "ksp_type": "preonly",
                # "ksp_type": "gmres",
                "pc_type": "lu",
                "pc_factor_mat_solver_type": "mumps",
                "mat_mumps_icntl_33": 0,
            },
            "newton_verbose": False,
            "ksp_verbose": False,
            "debug": False,
            "linear_solver": "mumps",
            # "linear_solver": "gmres",
            "error_on_nonconvergence": True,
            "relative_tolerance": 1e-5,
            "absolute_tolerance": 1e-5,
            "maximum_iterations": 20,
            "report": False,
            "krylov_solver": {
                "nonzero_initial_guess": True,
                "absolute_tolerance": 1e-10,
                "relative_tolerance": 1e-10,
                "maximum_iterations": 1000,
                "monitor_convergence": False,
            },
            "lu_solver": {"report": False, "symmetric": False, "verbose": False},
        }

    def converged(self, r, p, i):
        self._converged_called = True

        # breakpoint()
        res = r.norm("l2")
        print(f"Mechanics solver residual: {res}")
        # if res > 0.1:
        #     breakpoint()
        if i > 17:
            breakpoint()

        if self.debug:
            if not hasattr(self, "_datacollector"):
                print("No datacollector registered with the NewtonSolver")

            else:
                self._residuals.append(res)

        return super().converged(r, p, i)

    def solver_setup(self, A, J, p, i):
        self._solver_setup_called = True
        super().solver_setup(A, J, p, i)

    def solve(self, t0: float, dt: float):
        self.t0 = t0
        self.dt = dt
        print("Solving mechanics")
        self._solve_called = True
        self.ode.step(t0, dt)
        # self._diff.zero()
        # self._state.vector().zero()
        # self._diff.axpy(1.0, self._state.vector())
        # self._diff.axpy(-1.0, self._prev_state)
        # print(f"NORM : {self._diff.norm('l2')}")
        # if self._diff.norm("l2") < 1e-8:
        #     print("Converged")
        #     return 0
        ret = super().solve(self._problem, self._state.vector())
        # self._state.vector().apply("insert")
        # self._prev_state.zero()
        # self._prev_state.axpy(1.0, self._state.vector())
        print("Done solving mechanics")
        # self.save_residuals()
        return ret

    # DEBUGGING
    # This is just to check if we are using the overloaded functions
    def check_overloads_called(self):
        assert getattr(self, "_converged_called", False)
        assert getattr(self, "_solver_setup_called", False)
        assert getattr(self, "_update_solution_called", False)
        assert getattr(self, "_solve_called", False)

    def update_solution(self, x, dx, rp, p, i):
        self._update_solution_called = True

        # Update x from the dx obtained from linear solver (Newton iteration) :
        # x = -rp*dx (rp : relax param)
        print(f"Updating mechanics solution with relax parameter {rp}, iteration {i}")
        # self.odels
        # .step(t0=self.t0, dt=self.dt)
        # breakpoint()

        super().update_solution(x, dx, rp, p, i)
        if self._update_cb is not None:
            self._update_cb(x)


class MechanicsProblem(pulse.MechanicsProblem):
    def __init__(self, *args, **kwargs):
        self.ode = kwargs.pop("ode", None)
        super().__init__(*args, **kwargs)

    def _init_solver(self):
        if hasattr(self, "_dirichlet_bc"):
            bcs = self._dirichlet_bc
        else:
            bcs = []
        self._problem = NonlinearProblem(
            J=self._jacobian,
            F=self._virtual_work,
            bcs=bcs,
        )

        self.solver = NewtonSolver(
            problem=self._problem,
            state=self.state,
            # update_cb=self.material.active.update_prev,
            parameters=self.solver_parameters,
            ode=self.ode,
        )

    def solve(self, t0: float, dt: float):
        return self.solver.solve(t0, dt)
