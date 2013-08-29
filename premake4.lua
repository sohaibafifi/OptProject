-- A solution contains projects, and defines the available configurations
_ACTION = _ACTION or "gmake"

solution "OptProject"
   configurations { "Generic", "Debug", "Server" , "myMachine" }
   project "optproject"
        location "build"
        targetdir "build"
        kind "StaticLib"
        objdir "build/obj"
        language "C++"
        files { "**.h", "**.hpp", "**.c", "**.cpp" }
        links  {"m", "rt"}
        configuration "myMachine"
                defines { "NDEBUG"}
                flags { "OptimizeSpeed", "FloatFast" }
                buildoptions {"-Ofast", "-march=corei7", "-mtune=corei7", "-ffast-math", "-msse", "-msse2", "-msse3", "-pipe"}

        configuration "Generic"
                        defines { "NDEBUG"}
                        flags { "OptimizeSpeed", "FloatFast" }
                        buildoptions {"-Wall", "-Wno-sign-compare"}
                        linkoptions { "-static", "-static-libstdc++" }

        configuration "Server"
                defines { "NDEBUG"}
                flags { "OptimizeSpeed","StaticRuntime", "FloatFast" }
                buildoptions {"-w", "-march=core2", "-mtune=core2", "-pedantic", "-Wunused-parameter" }
                linkoptions { "-static", "-static-libstdc++" }

        configuration "Debug"
                defines {"DEBUG","_DEBUG", "TRACE"}
                flags { "Symbols" }
                buildoptions {"--coverage", "-fprofile-arcs", "-ftest-coverage" }
                links  {"gcov"}
