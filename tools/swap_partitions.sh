#!/bin/bash
find $PWD -name slurm\*sh -exec sed -i "s/sapphire/HighMem/g" {} \;
