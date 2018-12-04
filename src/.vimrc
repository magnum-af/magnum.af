" F2: Run Doxygen
nnoremap <F2> :!( cd .. && doxygen .doxygen-config && cd latex && make )<CR>
" F3: make clean 
nnoremap <F3> :! cd ../build  && make clean && cd - <cr>
" F4: cmake .. 
nnoremap <F4> :! cd ../build  && cmake .. && cd - <cr>
" F5: make
nnoremap <F5> :! cd ../build && make -j<cr>
" F6: run cpu
nnoremap <F6> :!time ../bin/magnum.af-cpu<cr>
" F7: run opencl
nnoremap <F7> :!time ../bin/magnum.af-opencl<cr>
" F9: Apply YCM FixIt
nnoremap <F9> :YcmCompleter FixIt<CR>
" F10: Goto definition with F3
nnoremap <F10> :YcmCompleter GoTo<CR>

" Using makeprg does not provide color output:
"nnoremap <F5> :make!<cr>
"set makeprg=make\ -C\ ../build\ -j
"set makeprg=[[\ -f\ Makefile\ ]]\ &&\ make\ \\\|\\\|\ make\ -C\ .. 
