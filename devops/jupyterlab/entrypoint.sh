#!/bin/bash
source /miniconda/etc/profile.d/conda.sh
conda activate waterkit

# Keep the container running in the foreground
tail -f /dev/null