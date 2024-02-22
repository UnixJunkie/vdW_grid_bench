#!/bin/bash

make

time ./grid_init -dx 0.1 -l data/lig.pqrs -p data/prot.pqrs
