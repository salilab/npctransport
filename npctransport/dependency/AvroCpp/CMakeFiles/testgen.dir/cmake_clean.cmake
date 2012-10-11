FILE(REMOVE_RECURSE
  "CMakeFiles/testgen"
  "testgen.hh"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/testgen.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
