" Cleanbuild and run all
nnoremap <F3> :! ./makeandrun.sh ../../.. <cr>
" Cleanbuild and run
nnoremap <F4> :! ./maketests.sh ../../.. && ./bin/%:r <cr>
" Build and run cpu
nnoremap <F5> :! cd build && make -j && cd .. && ./bin/test_cpu_%:r <cr>
" Build and run opencl
nnoremap <F6> :! cd build && make -j && cd .. && ./bin/test_opencl_%:r <cr>
" Build and run cuda
nnoremap <F7> :! cd build && make -j && cd .. && ./bin/test_cuda_%:r <cr>
