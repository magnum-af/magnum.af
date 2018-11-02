" build
nnoremap <F5> :make!<cr>
" run
nnoremap <F6> :!../bin/magnum.af-cpu<cr>

set makeprg=make\ -C\ ../build\ -j
"set makeprg=[[\ -f\ Makefile\ ]]\ &&\ make\ \\\|\\\|\ make\ -C\ .. 

" Apply YCM FixIt
nnoremap <F9> :YcmCompleter FixIt<CR>
" Goto definition with F3
nnoremap <F3> :YcmCompleter GoTo<CR>
