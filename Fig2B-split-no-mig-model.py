# Original script by JF McLaughlin; updated by CC Oliver Brown Mar 2024
# Updated for dadi v.2.1.1, Python 3.8.8
# The basis for this script is from Kevin Hawkins, 5 Feb 2016. I added in a few comments.
# Note that it differs from Lavretsky et al. 2015 Molecular Ecology, scaup, from Dryad.
# Numpy is the numerical library dadi is built upon
from numpy import array

import dadi

# import the demographic model. 
import Demographics2Dmod

# Load the data
data = dadi.Spectrum.from_file('nu1-vs-nu2-folded-spectrum.fs')
ns = data.sample_sizes
data.mask_corners()

# print number of samples to verify correct load
print(ns)

# Grid point settings will be used for extrapolation.
# Grid points need to be formated [n,n+10,n+20]. 
# Needs to be bigger than the number of samples you have (n>ns) and this will be a strong determination as to how long your program will run.

pts_l = [50,60,70]

# Call particular model to run, the model choosen here is split.w.migration
func = dadi.Demographics2D.split_mig

#    params = (nu1,nu2,T,m)
#    ns = (n1,n2)
#
#    Split into two populations of specifed size, with migration.
#    nu1: Size of population 1 after split.
#    nu2: Size of population 2 after split.
#    T: Time in the past of split (in units of 2*Na generations) 
#    m: Migration rate between populations (2*Na*m), but for this model m = 0
#    n1,n2: Sample sizes of resulting Spectrum
#    pts: Number of grid points to use in integration.

# Now let's optimize parameters for this model.

# The upper_bound and lower_bound lists are for use in optimization.
# Occasionally the optimizer will try wacky parameter values. We in particular
# want to exclude values with very long times, very small population sizes, or
# very high migration rates, as they will take a long time to evaluate.
# Parameters are: (nu1,nu2,T,m)

#Set the upper and lower bounds to make sure that the bounderies are 
#there. Suggested time parameters: lower 0, upper 5, migration 
#parameters: lower 0, upper 10,size parameters: lower 1e-2, upper 100 
upper_bound = [5, 5, 1, 0]
lower_bound = [1e-1, 1e-1, 0, 0]

# This is our initial guess for the parameters, which is somewhat arbitrary.
p0 = [1,1,1,1]
# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)

# Perturb our parameters before optimization. This does so by taking each
# parameter up to a factor of two up or down.
p0 = dadi.Misc.perturb_params(p0, fold=2, upper_bound=upper_bound,
                              lower_bound=lower_bound)

print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l, 
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=1)
# The verbose argument controls how often progress of the optimizer should be
# printed. It's useful to keep track of optimization process.
print('Finshed optimization **************************************************')

print('Best-fit parameters: {0}'.format(popt))

# Calculate the best-fit model AFS.
model = func_ex(popt, ns, pts_l)

# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, data)
print('Maximum log composite likelihood: {0}'.format(ll_model))

# The optimal value of theta given the model.
theta = dadi.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))
print(pts_l, upper_bound, lower_bound)

# Plot a comparison of the resulting fs with the data.
import pylab
pylab.figure(1)
dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=0.1, resid_range=1,
                                    pop_ids =('pop-nu1', 'pop-nu2's))

# This ensures that the figure pops up. It may be unecessary if you are using
# ipython.
pylab.show()
