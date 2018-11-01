" build
nnoremap <F5> :make!<cr>
" run
nnoremap <F6> :!../bin/magnum.af-cpu<cr>
" Apply YCM FixIt
map <F9> :YcmCompleter FixIt<CR>

set makeprg=make\ -C\ ../build\ -j
"set makeprg=[[\ -f\ Makefile\ ]]\ &&\ make\ \\\|\\\|\ make\ -C\ .. 
