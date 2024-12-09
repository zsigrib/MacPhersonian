#pragma once

// Describes how verbose a certain program (or complicated
// function) should be while it is running.
//
// - `silent` means that nothing should be printed,
// - `result` means only the result of the function should
//      be quickly summarized,
// - `info` means that the program should print status updates
//      about the main tasks it is doing, but not about progress
//      on these main tasks,
// - `checkpoints` means that the program should additionally
//      print a progress report on its main tasks, so it is clear
//      it is not stuck.
enum verboseness {silent, result, info, checkpoints};