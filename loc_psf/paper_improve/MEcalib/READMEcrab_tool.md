#Crab tools
有效面积计算工具路径：
/sharefs/hbkg/user/cwang/prepare_paper/tools_code/effect_area2.py
crab入射计算工具路径：
/sharefs/hbkg/user/cwang/prepare_paper/tools_code/calculate.py  crab_flux_absorb(energy_min,energy_max)

Energy-tau第一列是能量，第二列是tau，第100行是1keV和tau=2.4220002，也就是对于1keV的吸收是exp(-tau)

Energy-tau是对应柱密度为 1*10^22时的tau，对于6*10^21，tau = 6*10^21 / 10^22 * tau
考虑吸收，crab的谱为：9.8 * E ^-2.11 * exp(-tau)


Python插值：np.interp(lc_time,att_time1,att_ra)
