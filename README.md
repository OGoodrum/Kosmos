# Kosmos

## The **Kosmos-methoxy radical** problem

You are hired by an interstellar civilization of much greater intelligence that
ours (!!) who are seeding the universe with the building blocks needed for a
one of their major projects. The specific task you have agreed to help solve
for them is to manage the chemical reaction needed to form **methoxy radicals**
(which is made up of one carbon atom, one oxygen atom, and three hydrogen atoms).

Despite the civilization's greater intelligence, they are having trouble getting
the carbon, oxygen, and hydrogen atoms to combine correctly due to some serious synchronization
problems. They are able to create each atom (i.e., one atom equals one thread), and so 
the challenge is to get one carbon thread, three hydrogen threads, and one oxygen thread all together at the
same time, regardless of the number or order of thread/atom creation.

Each carbon atom invokes a function named `c_ready()` when it is ready to react;
each hydrogen atom invokes a function named `h_ready()` when it is ready to react; and finally each oxygen atom invokes a function named `o_ready()` when it is ready to react.
For this problem you are to complete the code for `c_ready`, `h_ready`, and `o_ready`, adding 
whatever code and data structures are needed to solve the problem.

Here are some other requirements for solutions:

* The `h_ready`, `c_ready`, and `o_ready` functions must only delay until there are at least
  one carbon, one oxygen, and three hydrogen atoms ready to combine into the radical.

* When an arriving atom/thread determines that a radical can be made,
then that atom/thread must cause the radical to be made to indicate
the identities of the carbon atom, oxygen atom, and hydrogen atoms.
Below is an example showing the format of the report, where each line
indicates the atoms in a radical, and the atom/thread in parentheses
was the one that initiated the creation of the radical. For example,
below the 7th radical consists of carbon atom 7, oxygen atom 7,
hydrogen atoms 21 & 22 & 23, and it was the carbon atom that triggered
the creation of the radical.

```
000: C000 O000 H000 H001 H002 (H000)
001: C001 O001 H003 H004 H005 (H003)
002: C002 O002 H006 H007 H008 (H006)
003: C003 O003 H009 H010 H011 (H010)
004: C004 O004 H012 H013 H014 (H013)
005: C005 O005 H015 H016 H017 (H017)
006: C006 O006 H018 H019 H020 (H019)
007: C007 O007 H021 H022 H023 (C007)
008: C008 O008 H024 H025 H026 (H024)
009: C009 O009 H027 H028 H029 (O009)
```

* Because there are several layers of scheduling involved with POSIX
threads and the Linux thread library, I am unable to reason
correctly about the fine details of program behavior by looking at
`printf` statements. Therefore all reports of radical creation are
actually written to an internal log, and that log is output at the end
of the program. You have been provided with the log routines (in
`logging.[ch]`), which include the use of a mutex to ensure there are
no race conditions when multiple radicals are formed at the same time.
