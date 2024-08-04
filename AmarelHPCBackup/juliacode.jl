using Distributed
using ProgressMeter
using KPM

@showprogress 1 "Computing..." for i in 1:50
    sleep(0.1)
end

