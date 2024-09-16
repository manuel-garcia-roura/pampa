#include "utils.hxx"

/* Remove a directory: */
void utils::remove_directory(const std::string& dir) {
   
   /* Remove the directory if it exists: */
   struct stat st = {0};
   if (stat(dir.c_str(), &st) == 0)
      remove(dir.c_str());
   
}

/* Create a directory: */
void utils::create_directory(const std::string& dir) {
   
   /* Create the directory if it doesn't exist: */
   struct stat st = {0};
   if (stat(dir.c_str(), &st) == -1)
      mkdir(dir.c_str(), 0700);
   
}
