#!/bin/bash
cd "$(dirname "${BASH_SOURCE[0]}")"
../../magnum.af dwp_a_jump_nonequi.py test_run_x x
../../magnum.af dwp_a_jump_nonequi.py test_run_y y
../../magnum.af dwp_a_jump_nonequi.py test_run_z z
