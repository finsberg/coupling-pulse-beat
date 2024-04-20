import dolfin
import pulse
import logging
import ufl_legacy as ufl
import numpy as np
from dataclasses import dataclass, field
from typing import Callable
import numpy.typing as npt
from collections import deque

logger = logging.getLogger(__name__)


def enlist(obj):
    try:
        return list(obj)
    except TypeError:
        return [obj]


@dataclass
class ODESolver:  # (beat.odesolver.DolfinODESolver):
    Ta: dolfin.Function
    states: npt.NDArray
    parameters: npt.NDArray
    fun: Callable
    fun_monitor: Callable
    Ta_index: int
    lmbda_index: int
    dLambda_index: int
    f0: ufl.core.expr.Expr
    missing_variables: npt.NDArray | None = None
    num_missing_variables: int = 0
    _kwargs: dict[str, npt.NDArray] = field(default_factory=dict)

    def __post_init__(self):
        if self.missing_variables is not None:
            self._kwargs["missing_variables"] = self.missing_variables
        self.next_states = self.states.copy()
        self.prev_Ta = dolfin.Function(self.Ta.function_space())
        self.lmbda = dolfin.Function(self.activation_space)
        self.prev_lmbda = dolfin.Function(self.activation_space)
        self.prev_prev_lmbda = dolfin.Function(self.activation_space)
        self.lmbda.vector()[:] = self.parameters[self.lmbda_index]

        self.dt = 0.0
        self.t = 0.0
        self.update_prev_Ta()
        self.update_prev_lmbda()

    @property
    def dLambda(self):
        if self.dt < 1e-8:
            return 0.0

        value = (
            self.lmbda.vector().get_local() - self.prev_lmbda.vector().get_local()
        ) / self.dt

        value = np.where(np.abs(value) < 1e-9, 0.0, value)
        return value

    def update_prev_Ta(self):
        self.prev_Ta.vector().zero()
        self.prev_Ta.vector().axpy(1.0, self.Ta.vector())

    def update_prev_lmbda(self):
        print("UPDATE PREV LAMBDA")
        self.prev_prev_lmbda.vector().zero()
        self.prev_prev_lmbda.vector().axpy(1.0, self.prev_lmbda.vector())
        self.prev_lmbda.vector().zero()
        self.prev_lmbda.vector().axpy(1.0, self.lmbda.vector())

    @property
    def activation_space(self) -> dolfin.FunctionSpace:
        return self.Ta.function_space()

    @property
    def num_points(self) -> int:
        return self.states.shape[1]

    @property
    def num_states(self) -> int:
        return self.states.shape[0]

    def step(self, t0: float, dt: float) -> None:
        # print("Old states: ", self.states[:, 0])
        self.next_states[:] = self.fun(
            states=self.states, t=t0, parameters=self.parameters, dt=dt, **self._kwargs
        )
        # print("New states: ", self.next_states[:, 0])
        # self.update_Ta(t0 + dt)

    @property
    def shape_missing_values(self) -> tuple[int, int]:
        return (self.num_missing_variables, self.num_points)

    @property
    def shape(self) -> tuple[int, int]:
        return (self.num_states, self.num_points)

    def update_Ta(self, t):
        self.t = t
        monitor = self.fun_monitor(
            t=t,
            states=self.next_states,
            parameters=self.parameters,
            missing_variables=self.missing_variables,
        )
        Ta = monitor[self.Ta_index]
        self.Ta.vector().set_local(np.where(Ta < 0, 0, Ta))

    def update_states(self):
        self.states[:] = self.next_states
        self.update_prev_Ta()
        self.update_prev_lmbda()

    def update_stretch(self, F, dt):
        self.dt = dt
        expr = ufl.sqrt(ufl.inner(F * self.f0, F * self.f0))
        new_lmbda = dolfin.project(expr, self.activation_space)

        self.lmbda.vector().zero()
        self.lmbda.vector().axpy(1.0, new_lmbda.vector())

        self.parameters[self.lmbda_index] = self.lmbda.vector().get_local()
        self.parameters[self.dLambda_index] = self.dLambda

    def reset(self):
        self.Ta.vector().zero()
        self.Ta.vector().axpy(1.0, self.prev_Ta.vector())
        self.next_states[:] = self.states
        self.lmbda.vector().zero()
        self.lmbda.vector().axpy(1.0, self.prev_lmbda.vector())
        self.prev_lmbda.vector().zero()
        self.prev_lmbda.vector().axpy(1.0, self.prev_prev_lmbda.vector())


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


class ContinuationProblem(NonlinearProblem):
    def __init__(self, J, F, bcs, **kwargs):
        # self.problem = problem
        super().__init__(J, F, bcs, **kwargs)
        self._J = J
        self._F = F

        self.bcs = [bcs]

        # super(ContinuationProblem, self).__init__()

        self.fres = deque(maxlen=2)

        self.first_call = True
        self.skipF = False

        self._assemble_jacobian = True

    def form(self, A, P, b, x):
        # pb = self.problem
        if self._assemble_jacobian:
            dolfin.assemble_system(self._J, self._F, self.bcs, A_tensor=A, b_tensor=b)
        else:
            dolfin.assemble(self._F, tensor=b)
            if self.bcs:
                for bc in self.bcs:
                    bc.apply(b)
        self._assemble_jacobian = not self._assemble_jacobian

        return
        # Timer("ContinuationSolver: form")
        # pb = self.problem

        # # check if we need to assemble the jacobian
        # if self.first_call:
        #     reset_jacobian = True
        #     self.first_call = False
        #     self.skipF = True
        # else:
        #     reset_jacobian = b.empty() and not A.empty()
        #     self.skipF = reset_jacobian

        #     if len(self.fres) == 2 and reset_jacobian:
        #         if self.fres[1] < 0.1 * self.fres[0]:
        #             debug("REUSE J")
        #             reset_jacobian = False

        # if reset_jacobian:
        #     # requested J, assemble both
        #     debug("ASSEMBLE J")
        #     assemble_system(pb.dG, pb.G, pb.bcs, x0=x, A_tensor=A, b_tensor=b)

    def J(self, A, x):
        pass

    def F(self, b, x):
        return
        # if self.skipF:
        #     return
        # pb = self.problem
        # assemble(pb.G, tensor=b)
        # for bc in pb.bcs:
        #     bc.apply(b)
        # self.fres.append(b.norm("l2"))


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
        self.dx = dolfin.Measure("dx", domain=state.function_space().mesh())
        self.volume = dolfin.assemble(dolfin.Constant(1) * self.dx)
        self.ode = ode
        self._problem = problem
        self._state = state
        self._update_cb = update_cb
        self._prev_state = dolfin.Vector(state.vector().copy())
        self._diff = dolfin.Vector(state.vector().copy())
        #
        # Initializing Newton solver (parent class)
        self.petsc_solver = dolfin.PETScKrylovSolver()
        super().__init__(
            self._state.function_space().mesh().mpi_comm(),
            self.petsc_solver,
            dolfin.PETScFactory.instance(),
        )
        self._handle_parameters(parameters)
        # self.dt_mech = 10.0
        # self.t_mech = 0.0
        # self.parameters["maximum_iterations"] = 50

    def _handle_parameters(self, parameters):
        # Setting default parameters
        params = type(self).default_solver_parameters()

        if parameters is not None:
            params.update(parameters)

        for k, v in params.items():
            if self.parameters.has_parameter(k):
                self.parameters[k] = v
            if self.parameters.has_parameter_set(k):
                for subk, subv in params[k].items():
                    self.parameters[k][subk] = subv
        petsc = params.pop("petsc", {})
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
        # self.parameters["relaxation_parameter"] = 0.8

    # def register_datacollector(self, datacollector):
    #     self._datacollector = datacollector

    @staticmethod
    def default_solver_parameters():
        return {
            "petsc": {
                "ksp_type": "preonly",
                # "ksp_type": "gmres",
                # "pc_type": "lu",
                "pc_type": "cholesky",
                "pc_factor_mat_solver_type": "mumps",
                "mat_mumps_icntl_33": 0,
                "mat_mumps_icntl_7": 6,
            },
            "newton_verbose": False,
            "ksp_verbose": False,
            "debug": False,
            "linear_solver": "gmres",
            # "preconditioner": "lu",
            # "linear_solver": "mumps",
            "error_on_nonconvergence": False,
            "relative_tolerance": 1e-5,
            "absolute_tolerance": 1e-5,
            "maximum_iterations": 100,
            "report": False,
            "krylov_solver": {
                "nonzero_initial_guess": False,
                "absolute_tolerance": 1e-13,
                "relative_tolerance": 1e-13,
                "maximum_iterations": 1000,
                "monitor_convergence": False,
            },
            "lu_solver": {"report": False, "symmetric": True, "verbose": False},
        }

    def converged(self, r, p, i):
        self._converged_called = True

        # breakpoint()
        res = r.norm("l2")
        print(f"Mechanics solver residual: {res}")
        # if res > 0.1:
        #     breakpoint()
        # if i > 17:
        #     breakpoint()

        if self.debug:
            if not hasattr(self, "_datacollector"):
                print("No datacollector registered with the NewtonSolver")

            else:
                self._residuals.append(res)

        return super().converged(r, p, i)

    def solver_setup(self, A, J, p, i):
        self._solver_setup_called = True
        super().solver_setup(A, J, p, i)

    def super_solve(self):
        return super().solve(self._problem, self._state.vector())

    def solve(self, t0: float, dt: float) -> tuple[int, bool]:
        self.t0 = t0
        self.dt = dt
        print("Solving mechanics")
        # print(
        #     "Missing variables: min: ",
        #     self.ode.missing_variables.min(),
        #     "max: ",
        #     self.ode.missing_variables.max(),
        # )
        self._solve_called = True

        # if self._diff.norm("l2") < 1e-8:
        #     print("Converged")
        #     return (0, True)
        # breakpoint()
        # print("ODE states min: ", self.ode.states.min(), "max: ", self.ode.states.max())
        # print(f"Before: ({t0}, {dt})", self._state.vector().get_local()[:10])

        # # self.update_stretch()
        # self.ode.step(self.t0, self.dt)
        # self.update_stretch()
        self.ode.step(self.t0, self.dt)
        self.ode.update_Ta(self.t0 + self.dt)
        try:
            nit, conv = self.super_solve()
        except:
            breakpoint()
        print("After: ", self._state.vector().get_local()[:10])
        # nit, conv = super().solve(self._problem, self._state.vector())
        if not conv:
            breakpoint()
        self.update_stretch()
        self.ode.update_states()

        self._diff.zero()
        # self._state.vector().zero()
        self._diff.axpy(1.0, self._state.vector())
        self._diff.axpy(-1.0, self._prev_state)
        print(f"Difference between previous state : {self._diff.norm('l2')}")

        # self.prev_lmbda.vector().zero()
        # self.prev_lmbda.vector().axpy(1.0, self.lmbda.vector())
        # self._state.vector().apply("insert")
        self._prev_state.zero()
        self._prev_state.axpy(1.0, self._state.vector())
        print("Done solving mechanics")
        # self.save_residuals()
        return (nit, conv)

    def reset(self):
        self._state.vector().zero()
        self._state.vector().axpy(1.0, self._prev_state)
        # self._prev_state.vector().zero()
        self.ode.reset()

    # DEBUGGING
    # This is just to check if we are using the overloaded functions
    def check_overloads_called(self):
        assert getattr(self, "_converged_called", False)
        assert getattr(self, "_solver_setup_called", False)
        assert getattr(self, "_update_solution_called", False)
        assert getattr(self, "_solve_called", False)

    @property
    def lmbda(self):
        return self.ode.lmbda

    def update_stretch(self):
        # Update stretch
        u = dolfin.split(self._state)[0]
        F = ufl.grad(u) + ufl.Identity(3)
        self.ode.update_stretch(F, self.dt)

    def update_solution(self, x, dx, rp, p, i):
        self._update_solution_called = True

        # Update x from the dx obtained from linear solver (Newton iteration) :
        # x = -rp*dx (rp : relax param)
        print(f"Updating mechanics solution with relax parameter {rp}, iteration {i}")

        super().update_solution(x, dx, rp, p, i)
        self.update_stretch()
        self.ode.step(self.t0, self.dt)
        self.ode.update_Ta(self.t0 + self.dt)

        if self._update_cb is not None:
            self._update_cb(x)


class MechanicsProblem(pulse.MechanicsProblem):
    def __init__(
        self,
        geometry,
        material,
        bcs=None,
        ode=None,
    ):
        logger.debug("Initialize mechanics problem")
        self.ode = ode
        self.geometry = geometry
        self.material = material
        self.bcs = bcs

        # # Make sure that the material has microstructure information
        # for attr in ("f0", "s0", "n0"):
        #     setattr(self.material, attr, getattr(self.geometry, attr))

        # self.solver_parameters = NonlinearSolver.default_solver_parameters()
        # if solver_parameters is not None:
        #     self.solver_parameters.update(**solver_parameters)

        self._init_spaces()
        self._init_forms()

    def _init_spaces(self):
        logger.debug("Initialize spaces for mechanics problem")
        mesh = self.geometry.mesh

        P2 = dolfin.VectorElement("Lagrange", mesh.ufl_cell(), 2)
        P1 = dolfin.FiniteElement("Lagrange", mesh.ufl_cell(), 1)

        # P2_space = FunctionSpace(mesh, P2)
        # P1_space = FunctionSpace(mesh, P1)
        self.state_space = dolfin.FunctionSpace(mesh, P2 * P1)

        self.state = dolfin.Function(self.state_space, name="state")
        self.state_test = dolfin.TestFunction(self.state_space)

    def _init_forms(self):
        logger.debug("Initialize forms mechanics problem")
        # Displacement and hydrostatic_pressure
        u, p = dolfin.split(self.state)
        v, q = dolfin.split(self.state_test)

        # Some mechanical quantities
        F = dolfin.variable(ufl.grad(u) + ufl.Identity(3))

        internal_energy = self.material.strain_energy(F) - p * (ufl.det(F) - 1)

        self._virtual_work = dolfin.derivative(
            internal_energy * self.geometry.dx,
            self.state,
            self.state_test,
        )
        # external_work = self._external_work(u, v)
        # if external_work is not None:
        #     self._virtual_work += external_work

        self._set_dirichlet_bc()

        self._jacobian = dolfin.derivative(
            self._virtual_work,
            self.state,
            dolfin.TrialFunction(self.state_space),
        )
        self._init_solver()

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
        self._problem = ContinuationProblem(
            J=self._jacobian,
            F=self._virtual_work,
            bcs=bcs,
        )

        self.solver = NewtonSolver(
            problem=self._problem,
            state=self.state,
            # update_cb=self.material.active.update_prev,
            # parameters=self.solver_parameters,
            ode=self.ode,
        )

    def solve(self, t0: float, dt: float):
        return self.solver.solve(t0, dt)
