" run cpu
nnoremap <F4> :! cd ../build && make -j && cd ../ <cr>
" run cpu
nnoremap <F6> :! ../bin/test_%:r_cpu <cr>
" run opencl
nnoremap <F7> :! ../bin/test_%:r_opencl <cr>
" run cuda
nnoremap <F8> :! ../bin/test_%:r_cuda <cr>
" run all
nmap <F5> <F4> <F6> <F7> <F8>
