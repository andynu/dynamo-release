.. automodule:: dynamo

API
===

Import dynamo as::

   import dynamo as dyn

Preprocessing (pp)
~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: _autosummary

   pp.cell_cycle_scores
   pp.recipe_monocle

Tools (tl)
~~~~~~~~~~

.. autosummary::
   :toctree: _autosummary

   tl.DDRTree
   tl.cell_velocities
   tl.cell_wise_confidence
   tl.cluster_field
   tl.confident_cell_velocities
   tl.converter
   tl.diffusion
   tl.dynamics
   tl.expected_return_time
   tl.gene_wise_confidence
   tl.generalized_diffusion_map
   tl.hdbscan
   tl.mnn
   tl.moments
   tl.neighbors
   tl.psl
   tl.reduceDimension
   tl.stationary_distribution


Estimation (est)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Conventional scRNA-seq (est.csc)**

.. autosummary::
   :toctree: _autosummary

   est.csc.fit_alpha_degradation
   est.csc.fit_alpha_synthesis
   est.csc.fit_f
   est.csc.fit_gamma_lsq
   est.csc.fit_linreg
   est.csc.sol_p
   est.csc.sol_s
   est.csc.sol_u
   est.csc.solve_f

**Time-resolved metabolic labeling based scRNA-seq (est.tsc)**

.. autosummary::
   :toctree: _autosummary

   est.tsc.Estimation_DeterministicKin
   est.tsc.Estimation_DeterministicKinNosp
   est.tsc.Estimation_MomentDegNosp
   est.tsc.Estimation_MomentKin
   est.tsc.Estimation_MomentKinNosp
   est.tsc.Mixture_KinDeg_NoSwitching
   est.tsc.kinetic_estimation

Vector field (vf)
~~~~~~~~~~~~~~~~~

.. autosummary:: 
   :toctree: _autosummary

   vf.Ao_pot_map
   vf.DiffusionMatrix
   vf.IntGrad
   vf.MFPT
   vf.Potential
   vf.SparseVFC
   vf.VectorField
   vf.Wang_LAP
   vf.Wang_action
   vf.acceleration
   vf.action
   vf.alignment
   vf.cell_accelerations
   vf.curl
   vf.curvature
   vf.diffusionMatrix
   vf.divergence
   vf.gen_fixed_points
   vf.gen_gradient
   vf.path_integral
   vf.rank_acceleration_genes
   vf.rank_curvature_genes
   vf.rank_divergence_genes
   vf.rank_speed_genes
   vf.solveQ
   vf.speed
   vf.topography
   vf.torsion
   vf.transition_rate
   vf.vfGraph

Prediction (pd)
~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: _autosummary

   pd.fate
   pd.fate_bias
   pd.state_graph

Plotting (pl)
~~~~~~~~~~~~~

.. autosummary::
   :toctree: _autosummary

   pl.basic_stats
   pl.cell_cycle_scores
   pl.cell_wise_vectors
   pl.curl
   pl.curvature
   pl.divergence
   pl.dynamics
   pl.exp_by_groups
   pl.fate_bias
   pl.feature_genes
   pl.grid_vectors
   pl.jacobian
   pl.jacobian_heatmap
   pl.jacobian_kinetics
   pl.kinetic_curves
   pl.kinetic_heatmap
   pl.line_integral_conv
   pl.nneighbors
   pl.pca
   pl.phase_portraits
   pl.plot_energy
   pl.plot_fixed_points
   pl.plot_flow_field
   pl.plot_nullclines
   pl.plot_separatrix
   pl.plot_traj
   pl.save_fig
   pl.scatters
   pl.show_fraction
   pl.show_landscape
   pl.speed
   pl.state_graph
   pl.streamline_plot
   pl.topography
   pl.trimap
   pl.tsne
   pl.umap
   pl.variance_explained

Moive (mv)
~~~~~~~~~~~~~

.. autosummary::
   :toctree: _autosummary

   mv.StreamFuncAnim
   mv.animate_fates

Simulation (sim)
~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: _autosummary

   sim.Gillespie
   sim.Simulator
   sim.Ying_model
   sim.evaluate
   sim.neurogenesis
   sim.state_space_sampler
   sim.toggle
   sim.two_genes_motif

External (ext)
~~~~~~~~~~~~~~

.. autosummary::
   :toctree: _autosummary

   ext.ddhodge
   ext.mutual_inform
   ext.scifate_glmnet
   ext.scribe

