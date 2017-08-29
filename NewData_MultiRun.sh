#!/bin/bash -l
#run it by: ./MultiRun.sh flilename(do not wirte dir here, write it in the FileDir)
#eg ./NewData_MultiRun.sh 'elg_240_sgc.v2.TSR.SSR.chunk21_subset' 'random-sweep.merged.chunk21_TSR_SSR_subset' 
export PYTHONPATH=.:$PYTHONPATH
Nproc=96    
RanFileNum=20 #number of sub_random files
#Njob=1+2*Nproc+Nproc*Nproc/2-Nproc/2    # 作业数目 DD+DR+RRself+RRcross
#FileDir="/n/des/kong.291/ELG_Data/ELG_New_Data/ELG_public_release_Data_Outputs/"
#SubFileDir="/n/des/kong.291/ELG_Data/ELG_New_Data/ELG_public_release_Data_Outputs/"
FileDir="/global/homes/h/huikong/eboss/ELG_Old_Data_Output/"
SubFileDir="/global/homes/h/huikong/eboss/ELG_Old_Data_Output/"
str_data=$FileDir$1
str_random=$SubFileDir$2
strfits=".fits"
#echo $Njob
export PYTHONPATH=.:$PYTHONPATH
#for((i=1; i<=$Nproc; i++)); do
#    str5=$SubFileDir${str%.fits}${str2}${i}${str3}
#    echo $str_random $str5
#   ./SprtFile $str_random $str5 $i $Nproc
#done
function PushQue {    
Que="$Que $1"
Nrun=$(($Nrun+1))
}
function GenQue {     
OldQue=$Que
Que=""; Nrun=0
for PID in $OldQue; do
if [[ -d /proc/$PID ]]; then
PushQue $PID
fi
done
}
function ChkQue {     
OldQue=$Que
for PID in $OldQue; do
if [[ ! -d /proc/$PID ]] ; then
GenQue; break
fi
done
}
i=0
./NewData_CorrFun $str_data${i}$strfits $str_data${i}$strfits 1 0 0 #DD
for((i=1; i<$RanFileNum; i++)); do
    ./NewData_CorrFun $str_data${i}$strfits $str_data${i}$strfits 1 $i $i & #DD
    PID=$!
    PushQue $PID
    while [[ $Nrun -ge $Nproc ]]; do
        ChkQue
        sleep 1
    done
done
for((i=0; i<$RanFileNum; i++)); do
    ./NewData_CorrFun $str_random${i}$strfits $str_random${i}$strfits 2 $i $i & #RR self
    PID=$!
    PushQue $PID
    while [[ $Nrun -ge $Nproc ]]; do
        ChkQue
        sleep 1
    done
done
for((i=0; i<$RanFileNum; i++)); do #DR
    for((j=0; j<$RanFileNum; j++)); do
        ./NewData_CorrFun $str_data${i}$strfits $str_random${j}$strfits 3 $i $j &
        PID=$!
        PushQue $PID
        while [[ $Nrun -ge $Nproc ]]; do
            ChkQue
            sleep 1
        done
    done
done
for((i=0; i<$RanFileNum; i++)); do #RR_cross
    for((j=i+1; j<$RanFileNum; j++)); do
            ./NewData_CorrFun $str_random${i}$strfits $str_random${j}$strfits 4 $i $j &
            PID=$!
            PushQue $PID
            while [[ $Nrun -ge $Nproc ]]; do
            ChkQue
            sleep 1
            done
     done
done
for((i=0; i<$RanFileNum; i++)); do #DD_cross
    for((j=i+1; j<$RanFileNum; j++)); do
            ./NewData_CorrFun $str_data${i}$strfits $str_data${j}$strfits 5 $i $j &
            PID=$!
            PushQue $PID
            while [[ $Nrun -ge $Nproc ]]; do
            ChkQue
            sleep 1
            done
     done
done
wait

exit
