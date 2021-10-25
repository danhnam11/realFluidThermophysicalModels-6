## 1. select the real-fluid models of interest (refer to available options)

## 2. run blockMesh to generate mesh:
blockMesh

## 2. Then run realFluidReactingFoam solver:
realFluidReactingFoam


# IMPORTANT NOTE:
You should using solution of 1step Mechanism as initial value instead of
running directly GRI 3.0 detail Mechanism from 0 to reduce the computational
time.

