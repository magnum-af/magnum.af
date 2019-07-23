" Run current file with python3
nnoremap <F5> :!python3 %<cr>

" Alternative causing ImportError: dynamic module does not define init function on some systems
nnoremap <F6> :exec '!python3' shellescape(@%, 1)<cr>
