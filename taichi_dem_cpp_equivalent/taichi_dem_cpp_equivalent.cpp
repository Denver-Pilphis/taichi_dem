// taichi_dem_cpp_equivalent.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
#include "pch.h"

//=======================================================================
// entrance
//=======================================================================
int main()
{
    DEMSolverConfig config;
    DEMSolver solver(config);
    solver.init_particle_fields("input.p4p");
    solver.init_grid_fields(grid_n);
    if (bond) solver.init_simulation();
    std::cout << "Grid size: " << grid_n << "x" << grid_n << std::endl;
    
    /*
    for (const std::filesystem::directory_entry& entry : std::filesystem::directory_iterator("./output_data"))
        std::filesystem::remove_all(entry.path());
    */
    try
    {
        std::filesystem::remove("output.p4p");
    }
    catch (std::exception ex)
    {}

    try
    {
        std::filesystem::remove("output.p4c");
    }
    catch (std::exception ex)
    {
    }

    Integer step = 0;
    for (; step < solver.config.nsteps; ++step)
    {
        std::cout << step << " ";
        solver.run_simulation();
        if (step % solver.config.saving_interval_steps == 0)
            solver.save("output", solver.config.dt * (Real)step);
    }
    // Final save
    solver.save("output", solver.config.dt * (Real)step);

    // Clear up all the datasets
    // Clear StructFieldObjs
    for (Integer i = 0; i < solver.gf.rows(); ++i)
        if (solver.gf[i])
        {
            delete solver.gf[i];
            solver.gf[i] = nullptr;
        }
    
    for (Integer i = 0; i < solver.wf.rows(); ++i)
        if (solver.wf[i])
        {
            delete solver.wf[i];
            solver.wf[i] = nullptr;
        }

    for (Integer i = 0; i < solver.cf.rows(); ++i)
        for (Integer j = 0; j < solver.cf.cols(); ++j)
        {
            if (i >= j) continue;
            if (solver.cf(i, j))
            {
                delete solver.cf(i, j);
                solver.cf(i, j) = nullptr;
            }
        }

    for (Integer i = 0; i < solver.wcf.rows(); ++i)
        for (Integer j = 0; j < solver.wcf.cols(); ++j)
        {
            if (i >= j) continue;
            if (solver.wcf(i, j))
            {
                delete solver.wcf(i, j);
                solver.wcf(i, j) = nullptr;
            }
        }

    return 0;
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
