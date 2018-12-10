" execute
nnoremap <F5> :! ../magnum.af main_infineon_CoFeB_stress.cpp $PWD/todel 0 ../vortex/llg/plot_hysteresis.sh <cr> 
" remove todel and execute into todel
nnoremap <F6> :!if [ -d "todel" ]; then rm -r todel/; fi;  ../magnum.af main_infineon_CoFeB_stress.cpp $PWD/todel/ 0 ../vortex/llg/plot_hysteresis.sh <cr> 
" remove todel 
nnoremap <F7> :!rm -r todel/ <cr>
