" cmake .. 
nnoremap <F4> :! cd ../build  && cmake .. && cd - <cr>
" make
nnoremap <F5> :make!<cr>
" run
nnoremap <F6> :!../bin/magnum.af-cpu<cr>

" run opencl
nnoremap <F7> :!../bin/magnum.af-opencl<cr>

set makeprg=make\ -C\ ../build\ -j
"set makeprg=[[\ -f\ Makefile\ ]]\ &&\ make\ \\\|\\\|\ make\ -C\ .. 

" Apply YCM FixIt
nnoremap <F9> :YcmCompleter FixIt<CR>
" Goto definition with F3
nnoremap <F3> :YcmCompleter GoTo<CR>
