import os

def write_sub_file(proc_id, queue_value):
    sub_file_content = f"""
Universe              = vanilla
executable            = Test.sh
arguments             = $(ProcId) {proc_id}
GetEnv                = True
output                = output/$(ClusterId).$(ProcId).out
error                 = error/$(ClusterId).$(ProcId).err
log                   = log/$(ClusterId).log
requirements          = (OpSysAndVer =?= "CentOS7")
+JobFlavour           = "tomorrow"

queue {queue_value}
"""
    sub_file_name = f"runTest_{proc_id}.sub"
    with open(sub_file_name, 'w') as sub_file:
        sub_file.write(sub_file_content)
    return sub_file_name

def submit_job(sub_file_name):
    os.system(f"condor_submit {sub_file_name}")

if __name__ == "__main__":
    # Define queue values for procID 0 and 1
    queue_value_0 = 2581
    queue_value_1 = 43

    sub_file_name_0 = write_sub_file(0, queue_value_0)
    submit_job(sub_file_name_0)

    sub_file_name_1 = write_sub_file(1, queue_value_1)
    submit_job(sub_file_name_1)
