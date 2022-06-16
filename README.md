# hxmt_scan

## 介绍
本仓库为HXMT扫描巡天数据的处理相关软件，\
仅限内部交流使用，暂未适应对外开放使用。\
当前包括ME巡天、HE巡天、PSF标定、PSF定位程序


## ME巡天使用说明
### 自动运行
2021年8月开始，在AFS账户luoqi节点hlogin06自动运行
```buildoutcfg
screen -S MEscan
source /sharefs/hbkg/user/luoqi/home/scan_MEgenlc_env.sh
python ./ME/ME_SCAN/genlc/multi_time_merun.py P0301240
```
退出：ctrl+A+D\
重连：screen -r



## HE巡天使用说明
### 自动运行
2021年8月开始，在AFS账户luoqi节点hlogin05自动运行
```buildoutcfg
screen -S HEscan
source /sharefs/hbkg/user/luoqi/home/mypython
python ./HE/HXMT_GPS_Promgram/HEsub_SCAN_Auto.py
```
退出：ctrl+A+D\
重连：screen -r

HE巡天的具体操作方法可参考./HE/HXMT_GPS_Promgram/readme.txt

## PSF标定使用说明

1.  xxxx
2.  xxxx
3.  xxxx


### 作业系统tips
根据实验需求需求，HTCondor 计算集群增加了几类特殊用途的作业队列，信息如下：

        1）长作业队列：用于运行超长时间的作业：
              - 单作业时长限制30天，集群长作业队队列的总运行作业数量为30个；
              - 提交作业示例：hep_sub job.sh –wt long

        2）短作业队列：用于运行较短时长的作业：
              - 单作业时长限制30分钟，集群短作业队列的总运行作业数量为40个
              - 用户短作业的优先级独立核算，即只按短作业的资源使用情况核算每个用户的短作业优先级
              - 提交示例：hep_sub job.sh –wt short

        3）超短作业队列：用于快速测试的作业：
              - 作业时长限制5分钟，集群超短作业队列的总运行作业数量为12个
              - 提交示例：hep_sub job.sh –wt test

#### 当前主要贡献人员
曩毅、赛娜、郭程程、廖进元、罗琦、李承奎


