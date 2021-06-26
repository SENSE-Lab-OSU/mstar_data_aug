function im_preproc=pre_process(img)
    im_preproc=squeeze(abs(img));
    im_preproc=im_preproc./norm(im_preproc,'fro');
    im_preproc=(20*log10(im_preproc));
end