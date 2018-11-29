" make clean 
nnoremap <F3> :! cd ../build  && make clean && cd - <cr>
" cmake .. 
nnoremap <F4> :! cd ../build  && cmake .. && cd - <cr>
" make
nnoremap <F5> :! cd ../build && make -j<cr>
"nnoremap <F5> :make!<cr>
"set makeprg=make\ -C\ ../build\ -j
"set makeprg=[[\ -f\ Makefile\ ]]\ &&\ make\ \\\|\\\|\ make\ -C\ .. 

" run
nnoremap <F6> :!../bin/magnum.af-cpu<cr>

" run opencl
nnoremap <F7> :!../bin/magnum.af-opencl<cr>

" Apply YCM FixIt
nnoremap <F9> :YcmCompleter FixIt<CR>
" Goto definition with F3
nnoremap <F10> :YcmCompleter GoTo<CR>
