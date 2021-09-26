1.程序和输出结果的迁移(可选):
原程序目录
/sharefs/hbkg/user/luoqi/GRB/work/ihep4/HXMT_GPS_Program

若需要迁移程序至dir1,输出结果目录放在dir2:
1)将HXMT_GPS_Program文件夹复制到dir1,并修改迁移后目录下的文件'dir_config.txt',
第一行改为dir1, 第二行改为dir2
2)修改完'dir_config.txt'后执行python mkdir_tree.py运行一次即可。
注：该步骤仅在每次修改目录时需要进行，日常运行时，不修改目录无需执行以上步骤。

2.端口执行程序:
产生光变和拟合当前分开，未集成在一起，因此执行分两部
0)source环境：
source /sharefs/hbkg/user/luoqi/home/mypython
如存在找不到headas等环境类型错误，可先执行以下语句再执行上一步：
source /sharefs/hbkg/user/luoqi/home/.bashrc
1)产生光变:
sh dir1/HXMT_GPS_Program/task_genlc.sh 观测号
示例:
sh /sharefs/hbkg/user/luoqi/GRB/work/ihep4/HXMT_GPS_Program/task_genlc.sh P010129506101
2)执行拟合:
sh dir1/HXMT_GPS_Program/task_lcfit.sh 观测号
示例:
sh /sharefs/hbkg/user/luoqi/GRB/work/ihep4/HXMT_GPS_Program/task_lcfit.sh P010129506101

3.提交作业执行程序(通常大批量观测号需处理时快速、方便):
产生光变和拟合当前分开，未集成在一起，因此执行分两部
0)source环境：
source /sharefs/hbkg/user/luoqi/home/mypython
如存在找不到headas等环境类型错误，可先执行以下语句再执行上一步：
source /sharefs/hbkg/user/luoqi/home/.bashrc
1)产生光变:
sh dir1/HXMT_GPS_Program/sub_gps_task.sh 观测号
示例:
sh /sharefs/hbkg/user/luoqi/GRB/work/ihep4/HXMT_GPS_Program/sub_gps_task.sh P010129506101
2)执行拟合:
***单个拟合时:
sh dir1/HXMT_GPS_Program/sub_gps_lcfit.sh 观测号
示例:
sh /sharefs/hbkg/user/luoqi/GRB/work/ihep4/HXMT_GPS_Program/sub_gps_lcfit.sh P010129506101
***批量光变已生成，需要拟合时:
sh dir1/gps_mission/lcfit_mission.sh
此文件记录了批量光变执行完的观测号，提交了对应拟合后请清空
示例:
sh /sharefs/hbkg/user/luoqi/GRB/work/ihep4/gps_mission/lcfit_mission.sh
rm /sharefs/hbkg/user/luoqi/GRB/work/ihep4/gps_mission/lcfit_mission.sh


####可能报错的环节####
0.以上所有程序基于python3环境,检查是否符合
1.文件授权问题，排除方法:对报错文件进行可读可写可执行授权
2.环境问题，排除方法：前述提供了一个可选环境方案，尝试后都不能解决时可自行完整安装heasoft，需要有grppha、xspec等程序包
3.执行作业时held, hep_q -u -hold查看，通常是内存超出导致，可以临时修改提高sub_gps_*.sh里-mem的数值后再次提交(表示兆字节内存)；或者端口执行

