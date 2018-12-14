" run python2 with argument todel/py1
nnoremap <buffer> <F3> :exec '!if [ -d "todel/py1/" ]; then rm -r todel/py1/; fi; PYTHONPATH=../../../../build/src python2' shellescape(@%, 1)' todel/py1/'<cr>
" run python3  with argument todel/py1
nnoremap <buffer> <F4> :exec '!if [ -d "todel/py1/" ]; then rm -r todel/py1/; fi; PYTHONPATH=../../../../build/src python3' shellescape(@%, 1)' todel/py1/'<cr>
" execute
nnoremap <F5> :! ../../../magnum.af main_infineon_CoFeB_stress.cpp $PWD/todel 0 plot.sh <cr> 
" remove todel and execute into todel
nnoremap <F6> :!if [ -d "todel" ]; then rm -r todel/; fi;  ../../../magnum.af main_infineon_CoFeB_stress.cpp $PWD/todel/ 0 plot.sh <cr> 
" remove todel 
nnoremap <F7> :!rm -r todel/ <cr>
