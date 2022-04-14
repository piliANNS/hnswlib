#include <iostream>
#include <fstream>
#include <queue>
#include <chrono>
#include "hnswlib/hnswlib.h"


#include <unordered_set>

using namespace std;
using namespace hnswlib;

// typedef unsigned char dataType;
// typedef int           resType;
typedef float dataType;
typedef float resType;

class StopW {
    std::chrono::steady_clock::time_point time_begin;
public:
    StopW() {
        time_begin = std::chrono::steady_clock::now();
    }

    float getElapsedTimeMicro() {
        std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
        return (std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_begin).count());
    }

    void reset() {
        time_begin = std::chrono::steady_clock::now();
    }

};



/*
* Author:  David Robert Nadeau
* Site:    http://NadeauSoftware.com/
* License: Creative Commons Attribution 3.0 Unported License
*          http://creativecommons.org/licenses/by/3.0/deed.en_US
*/

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))

#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif


/**
* Returns the peak (maximum so far) resident set size (physical
* memory use) measured in bytes, or zero if the value cannot be
* determined on this OS.
*/
static size_t getPeakRSS() {
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ((fd = open("/proc/self/psinfo", O_RDONLY)) == -1)
        return (size_t)0L;      /* Can't open? */
    if (read(fd, &psinfo, sizeof(psinfo)) != sizeof(psinfo))
    {
        close(fd);
        return (size_t)0L;      /* Can't read? */
    }
    close(fd);
    return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t) (rusage.ru_maxrss * 1024L);
#endif

#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L;          /* Unsupported. */
#endif
}


/**
* Returns the current resident set size (physical memory use) measured
* in bytes, or zero if the value cannot be determined on this OS.
*/
static size_t getCurrentRSS() {
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
        (task_info_t)&info, &infoCount) != KERN_SUCCESS)
        return (size_t)0L;      /* Can't access? */
    return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    /* Linux ---------------------------------------------------- */
    long rss = 0L;
    FILE *fp = NULL;
    if ((fp = fopen("/proc/self/statm", "r")) == NULL)
        return (size_t) 0L;      /* Can't open? */
    if (fscanf(fp, "%*s%ld", &rss) != 1) {
        fclose(fp);
        return (size_t) 0L;      /* Can't read? */
    }
    fclose(fp);
    return (size_t) rss * (size_t) sysconf(_SC_PAGESIZE);

#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return (size_t)0L;          /* Unsupported. */
#endif
}


static void
get_gt(unsigned *massQA, dataType *massQ, size_t vecsize, size_t qsize,
       size_t vecdim, vector<std::priority_queue<std::pair<resType, labeltype >>> &answers, size_t k) {


    (vector<std::priority_queue<std::pair<resType, labeltype >>>(qsize)).swap(answers);
    //DISTFUNC<int> fstdistfunc_ = l2space.get_dist_func();
    cout << qsize << "\n";
    for (int i = 0; i < qsize; i++) {
        for (int j = 0; j < k; j++) {
            answers[i].emplace(0.0f, massQA[100 * i + j]);
        }
    }
}

static float
test_approx(dataType *massQ, size_t vecsize, size_t qsize, HierarchicalNSW<resType> &appr_alg, size_t vecdim,
            vector<std::priority_queue<std::pair<resType, labeltype >>> &answers, size_t k) {
    size_t correct = 0;
    size_t total = 0;
    //uncomment to test in parallel mode:
    //#pragma omp parallel for
    for (int i = 0; i < qsize; i++) {

        std::priority_queue<std::pair<resType, labeltype >> result = appr_alg.searchKnn(massQ + vecdim * i, k);
        std::priority_queue<std::pair<resType, labeltype >> gt(answers[i]);
        unordered_set<labeltype> g;
        total += gt.size();

        while (gt.size()) {


            g.insert(gt.top().second);
            gt.pop();
        }

        while (result.size()) {
            if (g.find(result.top().second) != g.end()) {

                correct++;
            } else {
            }
            result.pop();
        }

    }
    return 1.0f * correct / total;
}

static void
test_vs_recall(dataType *massQ, size_t vecsize, size_t qsize, HierarchicalNSW<resType> &appr_alg, size_t vecdim,
               vector<std::priority_queue<std::pair<resType, labeltype >>> &answers, size_t k) {
    vector<size_t> efs;// = { 10,10,10,10,10 };
    for (int i = k; i < 30; i++) {
        efs.push_back(i);
    }
    for (int i = 30; i < 100; i += 10) {
        efs.push_back(i);
    }
    for (int i = 100; i < 500; i += 40) {
        efs.push_back(i);
    }
    for (size_t ef : efs) {
        appr_alg.setEf(ef);
        StopW stopw = StopW();

        float recall = test_approx(massQ, vecsize, qsize, appr_alg, vecdim, answers, k);
        float time_us_per_query = stopw.getElapsedTimeMicro() / qsize;

        cout << ef << "\t" << recall << "\t" << time_us_per_query << " us\n";
        if (recall > 1.0) {
            cout << recall << "\t" << time_us_per_query << " us\n";
            break;
        }
    }
}

inline bool exists_test(const std::string &name) {
    ifstream f(name.c_str());
    return f.good();
}

// load file. store format: (uint32_t)num, (uint32_t)dim, (data_T)num * dim.
template<typename data_T>
void LoadBinToArray(std::string& file_path, data_T *data_m, uint32_t nums, uint32_t dims, bool non_header = false){
    std::ifstream file_reader(file_path.c_str(), ios::binary);
    if (!non_header){
        uint32_t nums_r, dims_r;
        file_reader.read((char *) &nums_r, sizeof(uint32_t));
        file_reader.read((char *) &dims_r, sizeof(uint32_t));
        if ((nums != nums_r) || (dims != dims_r)){
            printf("Error, file %s is error, nums_r: %u, dims_r: %u\n", file_path.c_str(), nums_r, dims_r);
            exit(1);
        }
    }

    file_reader.read((char *) data_m, nums * dims * sizeof(data_T));
    file_reader.close();
    printf("Load %u * %u Data from %s done.\n", nums, dims, file_path.c_str());
}

void sift_test1B() {
	
	
	int subset_size_milllions = 10;
	int efConstruction = 40;
	int M = 16;
	

    size_t vecsize = subset_size_milllions * 1000000;

    size_t qsize = 10000;
    //size_t vecdim = 128;
    size_t vecdim = 96;
    char path_index[1024];
    // string path_q = "/home/nfs_data/hujingbo99/sift/query.public.10K.u8bin";
    // string path_data = "/home/nfs_data/hujingbo99/sift/sift10m/base.10m.u8bin";
    // sprintf(path_index, "sift1b_%dm_ef_%d_M_%d.bin", subset_size_milllions, efConstruction, M);
    // string path_gt = "/home/nfs_data/hujingbo99/sift/sift10m/groundtruth." + to_string(subset_size_milllions) + "m.bin";
     
    string path_q = "/home/nfs_data/hujingbo99/deep/query.public.10K.fbin";
    string path_data = "/home/nfs_data/hujingbo99/deep/deep10m/base.10m.fbin";
    sprintf(path_index, "deep1b_%dm_ef_%d_M_%d.bin", subset_size_milllions, efConstruction, M);
    string path_gt = "/home/nfs_data/hujingbo99/deep/deep10m/groundtruth." + to_string(subset_size_milllions) + "m.bin";

    // unsigned char *massb = new unsigned char[vecdim];
    
    int gt_maxnum = 100;

    cout << "Loading GT:\n";
    unsigned *massQA = new unsigned[qsize * gt_maxnum];
    LoadBinToArray<unsigned>(path_gt, massQA, qsize, gt_maxnum);
    // ifstream inputGT(path_gt, ios::binary);
    // 
    // for (int i = 0; i < qsize; i++) {
    //     int t;
    //     inputGT.read((char *) &t, 4);
    //     inputGT.read((char *) (massQA + 1000 * i), t * 4);
    //     if (t != 1000) {
    //         cout << "err";
    //         return;
    //     }
    // }
    // inputGT.close();
	
    cout << "Loading queries:\n";
    dataType *massQ = new dataType[qsize * vecdim];
    LoadBinToArray<dataType>(path_q, massQ, qsize, vecdim);
    // ifstream inputQ(path_q, ios::binary);

    // for (int i = 0; i < qsize; i++) {
    //     int in = 0;
    //     inputQ.read((char *) &in, 4);
    //     if (in != 128) {
    //         cout << "file error";
    //         exit(1);
    //     }
    //     inputQ.read((char *) massb, in);
    //     for (int j = 0; j < vecdim; j++) {
    //         massQ[i * vecdim + j] = massb[j];
    //     }

    // }
    // inputQ.close();

    dataType *massB = new dataType[vecsize * vecdim];
    LoadBinToArray<dataType>(path_data, massB, vecsize, vecdim);

    // unsigned char *mass = new unsigned char[vecdim];
    // ifstream input(path_data, ios::binary);
    int in = 0;
    //L2SpaceI l2space(vecdim);
    L2Space l2space(vecdim);


    HierarchicalNSW<resType> *appr_alg;
    if (exists_test(path_index)) {
        cout << "Loading index from " << path_index << ":\n";
        appr_alg = new HierarchicalNSW<resType>(&l2space, path_index, false);
        cout << "Actual memory usage: " << getCurrentRSS() / 1000000 << " Mb \n";
    } else {
        cout << "Building index:\n";
        appr_alg = new HierarchicalNSW<resType>(&l2space, vecsize, M, efConstruction);


        // input.read((char *) &in, 4);
        // if (in != 128) {
        //     cout << "file error";
        //     exit(1);
        // }
        // input.read((char *) massb, in);

        // for (int j = 0; j < vecdim; j++) {
        //     mass[j] = massb[j] * (1.0f);
        // }
        // mass = massB;
        appr_alg->addPoint((void *) (massB), (size_t) 0);
        int j1 = 0;
        StopW stopw = StopW();
        StopW stopw_full = StopW();
        size_t report_every = 100000;
#pragma omp parallel for
        for (int i = 1; i < vecsize; i++) {
            // unsigned char mass[128];
            int j2=0;
#pragma omp critical
            {

                // input.read((char *) &in, 4);
                // if (in != 128) {
                //     cout << "file error";
                //     exit(1);
                // }
                // input.read((char *) massb, in);
                // for (int j = 0; j < vecdim; j++) {
                //     mass[j] = massb[j];
                // }
                j1++;
                j2=j1;
                if (j1 % report_every == 0) {
                    cout << j1 / (0.01 * vecsize) << " %, "
                         << report_every / (1000.0 * 1e-6 * stopw.getElapsedTimeMicro()) << " kips " << " Mem: "
                         << getCurrentRSS() / 1000000 << " Mb \n";
                    stopw.reset();
                }
            }
            appr_alg->addPoint((void *) (massB + j2 * vecdim), (size_t) j2);


        }
        // input.close();
        cout << "Build time:" << 1e-6 * stopw_full.getElapsedTimeMicro() << "  seconds\n";
        appr_alg->saveIndex(path_index);
    }


    vector<std::priority_queue<std::pair<resType, labeltype >>> answers;
    size_t k = 10;
    cout << "Parsing gt:\n";
    get_gt(massQA, massQ, vecsize, qsize, vecdim, answers, k);
    cout << "Loaded gt\n";
    for (int i = 0; i < 1; i++)
        test_vs_recall(massQ, vecsize, qsize, *appr_alg, vecdim, answers, k);
    cout << "Actual memory usage: " << getCurrentRSS() / 1000000 << " Mb \n";
    return;


}
