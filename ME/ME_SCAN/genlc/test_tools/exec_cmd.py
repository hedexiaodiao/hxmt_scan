import subprocess###
def exec_cmd(cmd,my_env):
    """Run shell command"""
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,
                         env=my_env,
                         shell=True)

    log_content = p.communicate()[0]
    #print(log_content.decode('utf-8'))

    return p.returncode, log_content
