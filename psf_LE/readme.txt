1.假设如LE提供已合并好的光变：
当前目录执行autorun_psf.py 仪器名称(LE/ME/HE)，提交作业命令将写在以下文件中，
'./command_'+instru.lower()+'_'+time.strftime("%y%m%d")+'.sh'
举例，今天21年12月19日，执行：
python autorun_psf.py LE
将得到当前目录下./command_le_211219.sh

然后在bash执行：
source /sharefs/hbkg/user/luoqi/home/mypython
sh ./command_le_211219.sh
将自动提交作业进行拟合。
如果不想提交作业，将command_le_211219.sh中的
sub_task_psf.sh改成task_psf.sh，逐个命令在终端运行(同时运行多个将占用大量资源)

2.合并观测号文件方式：
如果文件命名方式改了，
Fitstogether2.py不需要更改，
但是Merge_obs.py需要修改
041crab.all填入要合并的观测号
需要根据路径文件名修改config_he.xml
