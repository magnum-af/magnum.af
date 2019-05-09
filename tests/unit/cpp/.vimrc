" Cleanbuild and run all
nnoremap <F3> :! ./makeandrun.sh ../../.. <cr>
" Cleanbuild and run
nnoremap <F4> :! ./maketests.sh ../../.. && ./bin/%:r <cr>
" Build and run
nnoremap <F5> :! cd build && make -j && cd .. && ./bin/%:r <cr>
