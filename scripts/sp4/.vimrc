nnoremap <buffer> <F4> :exec '!PYTHONPATH=../../build/python python3' shellescape(@%, 1)<cr>
nnoremap <buffer> <F5> :!../magnum.af main_sp4.cpp $PWD/todel 0 plot_sp4.sh<cr>
