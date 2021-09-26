echo "This is a new HXMTDAS version (V2)"
echo "Update: 2018 11 18, by Zhao Haisheng"
echo "herspgen,merspgen and lerspgen can be used to generate response file"
echo "CALDB has been updated"
echo "lebkgmap update: smooth background spectrum"

export HEADAS=/home/hxmt/hxmtsoft2/soft/install/x86_64-unknown-linux-gnu-libc2.12
source $HEADAS/headas-init.sh
#export HEADAS=/home/hxmt/zhaohaish/soft/install/x86_64-unknown-linux-gnu-libc2.12/
#source $HEADAS/headas-init.sh
export CALDBALIAS=/home/hxmt/hxmtsoft2/CALDB/software/tools/alias_config.fits
export CALDB=/home/hxmt/hxmtsoft2/CALDB/
export CALDBCONFIG=/home/hxmt/hxmtsoft2/CALDB/caldb.config

echo "pay attention:"
echo "HEADAS:"$HEADAS
echo "CALDB:"$CALDB

#export PATH="/hxmt/soft/Develop/anaconda2/bin:$PATH"
export REFPATH="/home/hxmt/gemy/work/HDPC/hxmtbkg/hxmt_bkgest_v2-0.7/refdata/"

echo "New backgoud process for test"
alias hebkgmap="python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/hebkgmap.py"
alias tmphebkgmap="python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/tmphebkgmap.py"
alias tmphebkgmap2="python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/tmphebkgmap_20190507.py"
alias tmphebkgmap_nocor="python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/tmphebkgmap_nocor.py"
alias mebkgmap="python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/mebkgmap.py"
alias lebkgmap="python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/lebkgmap.py"
alias tmplebkgmap_noehk="python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/tmplebkgmap_noehk.py"
echo "hebkgmap -h to find the usage"
echo "mebkgmap -h to find the usage"
echo "lebkgmap -h to find the usage"

echo "New gti selection for ME & LE"
echo "Update: 2018 10 19"
echo "legti le_recon.fits meoldgti.fits menewgti.fits"
echo "megti me_grade.fits meoldgti.fits menewgti.fits"
alias megti='python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/megti.py'
alias legti='python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/legti.py'
alias legti2='python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/legti_par2.py'
alias hegti='python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/hegti.py'
alias hprint_detid='python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/hprint_detid.py'
alias hphase_cal='python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/hphase_cal.py'
alias hen2pi='python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/hen2pi.py'

alias hhe_spec2pi="python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/hhe_spec2pi.py"

echo "The background of Insight-HXMT has been updated in 2019-05-31."
echo "The usage for lebkgmap, mebgkmap and hebkgmap does not changing."
echo "The updated version could be performed by source /hxmt/home/hxmtsoft2/hxmtsoft_v2.01.sh or source /hxmt/home/hxmtsoft2/hxmtsoft_v2.01.sh"
