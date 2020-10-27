#!/bin/bash

LBPM_INSTALL_DIR=../../tests

mpirun --oversubscribe -np 8 $LBPM_INSTALL_DIR/lbpm_random_force_simulator input.db
