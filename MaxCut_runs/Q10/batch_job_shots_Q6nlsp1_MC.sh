#! /bin/bash
#SBATCH --job-name=iZNEQ6N1c
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=2-00:00:00
#SBATCH --no-requeue
#SBATCH --output=out
#SBATCH --error=err
#SBATCH --qos=long
# #SBATCH -p shared
#SBATCH --nodelist=cn450
set -x                          # Output commands
set -e                          # Abort o-14n errors
echo "Checking:"
pwd
hostname
date
env | sort < /dev/null > outfiles/ENVIRONMENT
echo "Starting:"
module load matlab
source activate base

nohup taskset -c 0 matlab -r "shot_noise_budget('6', '0', '1', '1', '', '')" < /dev/null > outfiles/o1.out0 2> outfiles/o1.err0&
nohup taskset -c 1 matlab -r "shot_noise_budget('6', '0', '2', '1', '', '')" < /dev/null > outfiles/o2.out0 2> outfiles/o2.err0&
nohup taskset -c 2 matlab -r "shot_noise_budget('6', '0', '3', '1', '', '')" < /dev/null > outfiles/o3.out0 2> outfiles/o3.err0&
nohup taskset -c 3 matlab -r "shot_noise_budget('6', '0', '4', '1', '', '')" < /dev/null > outfiles/o4.out0 2> outfiles/o4.err0&
nohup taskset -c 4 matlab -r "shot_noise_budget('6', '0', '5', '1', '', '')" < /dev/null > outfiles/o5.out0 2> outfiles/o5.err0&
nohup taskset -c 5 matlab -r "shot_noise_budget('6', '0', '6', '1', '', '')" < /dev/null > outfiles/o6.out0 2> outfiles/o6.err0&
nohup taskset -c 6 matlab -r "shot_noise_budget('6', '0', '7', '1', '', '')" < /dev/null > outfiles/o7.out0 2> outfiles/o7.err0&
nohup taskset -c 7 matlab -r "shot_noise_budget('6', '0', '8', '1', '', '')" < /dev/null > outfiles/o8.out0 2> outfiles/o8.err0&
nohup taskset -c 8 matlab -r "shot_noise_budget('6', '0', '9', '1', '', '')" < /dev/null > outfiles/o9.out0 2> outfiles/o9.err0&
nohup taskset -c 9 matlab -r "shot_noise_budget('6', '0', '10', '1', '', '')" < /dev/null > outfiles/o10.out0 2> outfiles/o10.err0&
nohup taskset -c 10 matlab -r "shot_noise_budget('6', '0', '11', '1', '', '')" < /dev/null > outfiles/o11.out0 2> outfiles/o11.err0&
nohup taskset -c 11 matlab -r "shot_noise_budget('6', '0', '12', '1', '', '')" < /dev/null > outfiles/o12.out0 2> outfiles/o12.err0&
nohup taskset -c 12 matlab -r "shot_noise_budget('6', '0', '13', '1', '', '')" < /dev/null > outfiles/o13.out0 2> outfiles/o13.err0&
nohup taskset -c 13 matlab -r "shot_noise_budget('6', '0', '14', '1', '', '')" < /dev/null > outfiles/o14.out0 2> outfiles/o14.err0&
nohup taskset -c 14 matlab -r "shot_noise_budget('6', '0', '15', '1', '', '')" < /dev/null > outfiles/o15.out0 2> outfiles/o15.err0&
nohup taskset -c 15 matlab -r "shot_noise_budget('6', '0', '16', '1', '', '')" < /dev/null > outfiles/o16.out0 2> outfiles/o16.err0&
nohup taskset -c 16 matlab -r "shot_noise_budget('6', '0', '17', '1', '', '')" < /dev/null > outfiles/o17.out0 2> outfiles/o17.err0&
nohup taskset -c 17 matlab -r "shot_noise_budget('6', '0', '18', '1', '', '')" < /dev/null > outfiles/o18.out0 2> outfiles/o18.err0&
nohup taskset -c 18 matlab -r "shot_noise_budget('6', '0', '19', '1', '', '')" < /dev/null > outfiles/o19.out0 2> outfiles/o19.err0&
nohup taskset -c 19 matlab -r "shot_noise_budget('6', '0', '20', '1', '', '')" < /dev/null > outfiles/o20.out0 2> outfiles/o20.err0&
nohup taskset -c 20 matlab -r "shot_noise_budget('6', '0', '21', '1', '', '')" < /dev/null > outfiles/o21.out0 2> outfiles/o21.err0&
nohup taskset -c 21 matlab -r "shot_noise_budget('6', '0', '22', '1', '', '')" < /dev/null > outfiles/o22.out0 2> outfiles/o22.err0&
nohup taskset -c 22 matlab -r "shot_noise_budget('6', '0', '23', '1', '', '')" < /dev/null > outfiles/o23.out0 2> outfiles/o23.err0&
nohup taskset -c 23 matlab -r "shot_noise_budget('6', '0', '24', '1', '', '')" < /dev/null > outfiles/o24.out0 2> outfiles/o24.err0&
nohup taskset -c 24 matlab -r "shot_noise_budget('6', '0', '25', '1', '', '')" < /dev/null > outfiles/o25.out0 2> outfiles/o25.err0&
nohup taskset -c 25 matlab -r "shot_noise_budget('6', '0', '26', '1', '', '')" < /dev/null > outfiles/o26.out0 2> outfiles/o26.err0&
nohup taskset -c 26 matlab -r "shot_noise_budget('6', '0', '27', '1', '', '')" < /dev/null > outfiles/o27.out0 2> outfiles/o27.err0&
nohup taskset -c 27 matlab -r "shot_noise_budget('6', '0', '28', '1', '', '')" < /dev/null > outfiles/o28.out0 2> outfiles/o28.err0&
nohup taskset -c 28 matlab -r "shot_noise_budget('6', '0', '29', '1', '', '')" < /dev/null > outfiles/o29.out0 2> outfiles/o29.err0&
nohup taskset -c 29 matlab -r "shot_noise_budget('6', '0', '30', '1', '', '')" < /dev/null > outfiles/o30.out0 2> outfiles/o30.err0&

wait
echo "Stopping:"
date
