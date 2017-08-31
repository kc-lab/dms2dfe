def test_datasets():
    import subprocess
    import os
    from os.path import exists
    # from dms2dfe import pipeline
    repo_n='ms_datasets'
    repo_v='0.0.4'
    datasets=[
              # 'Firnberg_et_al_2014',
              'Olson_et_al_2014',
              # 'Melnikov_et_al_2014',
              ]
    datasets_dh='%s-%s' % (repo_n,repo_v)
    if not exists(datasets_dh):
        subprocess.call('wget https://github.com/rraadd88/%s/archive/v%s.zip' % (repo_n,repo_v),shell=True)
        subprocess.call('unzip v%s.zip' % (repo_v),shell=True)

    os.chdir('%s-%s/analysis' % (repo_n,repo_v))
    for prj_dh in datasets:
        print prj_dh
        if exists(prj_dh):
            subprocess.call('dms2dfe %s' % (prj_dh),shell=True)
        # break


