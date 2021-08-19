import jax
import jax.numpy as jnp

@jax.jit
def mask_cat(cat, mask):
    ra, dec = cat[0], cat[1]



