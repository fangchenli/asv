Understanding ASV step detection from the beginning
===================================================

ASV helps answer a practical question: after someone changes a program, does
the program do the same work faster or slower?

This guide starts with what is being measured. It then explains how ASV
interprets the measurements. You do not need to know statistics, programming,
or the Potts model. The arithmetic uses addition, subtraction, multiplication,
and division, and each new notation is explained before it is used.

.. contents:: On this page
   :local:
   :depth: 1

What a benchmark is
-------------------

A program does work. For example, it might put a list of one thousand names
in alphabetical order. We can ask a computer to do that task and measure how
long it takes.

A **benchmark** is a repeatable task used to measure a program's performance.
For a useful comparison, we keep the task and its input the same: sort the
same list, rather than compare sorting a short list with sorting a long one.

In this guide, performance means the time needed to finish the task. A smaller
time means faster performance. Benchmarks can measure other things, such as
memory use, but we will use time throughout the example.

We measure time in **milliseconds**. One millisecond is one thousandth of a
second. A task taking 10 milliseconds takes 0.010 seconds. The small numbers
are convenient for our example; the reasoning would also apply to a task
taking several seconds.

Why we measure several versions of a program
--------------------------------------------

People change programs to add features or fix problems. Each saved version
can behave a little differently. In a version-control system such as Git,
a saved change is called a **commit**. Think of commits as points in the
program's history.

ASV, short for airspeed velocity, is a tool that can run benchmarks against
different versions of a program and record the results. This lets us compare
how long the same task takes before and after changes to the program.

If a newer version takes longer to do the same work, that is a **performance
regression**, meaning a slowdown. If it takes less time, that is a speedup.
Finding a slowdown tells us something changed; it does not by itself tell
us which line of code caused it or whether the change was intentional.

Why the measurements are not perfectly steady
---------------------------------------------

Even without changing the program, repeated timings can differ. The computer
may be doing other work, or the conditions inside the computer may vary.
We might measure 10 milliseconds once, 10.1 next time, and 9.9 another time.

These small variations are often called **noise**. Here noise means variation
that we do not want to interpret as a lasting performance change. It is not
sound. A particular variation could still have a real cause; calling it noise
means we are treating it as background variation for this analysis.

ASV can repeat measurements and summarize them into a benchmark result.
For this guide, imagine that we already have one recorded timing for each
program version. The step detector receives that history of results. It is
not starting with the program's source code or every individual timing repeat.

The difficulty is deciding whether a higher number is ordinary variation
or evidence that the program has become consistently slower.

The measurements we will use
----------------------------

Suppose we have eight successive versions of a program. We number them
1 through 8 just to make the example easy to follow. These are example labels,
not real Git identifiers.

The recorded timings are::

    Version          1      2      3      4      5      6      7      8
    Measured time   10     10.1    9.9   10     12     12.1   11.9   12

All the times in that row are milliseconds. For example, version 2 took
10.1 milliseconds and version 6 took 12.1 milliseconds.

The first four numbers are close to 10. The last four are close to 12.
A reasonable explanation is that the task used to take about 10 milliseconds
and began taking about 12 milliseconds at version 5.

That explanation is something to infer from the measurements. We have not
told ASV in advance that the change happens at version 5.

What fitted means
-----------------

**Measured time** answers: what did we record when we timed the task?

**Fitted time** answers: what time does our proposed explanation assign to
this version? The analysis chooses these numbers to summarize the pattern
in the measurements. Fitted means chosen to match the data according to a
rule, which we will build up below.

For our example, one possible description is:

* Versions 1 through 4 take about 10 milliseconds.
* Versions 5 through 8 take about 12 milliseconds.

Writing that description underneath the measurements gives::

    Version          1      2      3      4      5      6      7      8
    Measured time   10     10.1    9.9   10     12     12.1   11.9   12
    Fitted time     10     10     10     10     12     12     12     12

At version 2, the measured time is 10.1 and the fitted time is 10. The fitted
row expresses the explanation that the extra 0.1 is ordinary variation around
a usual time of 10.

The fitted row is not another set of stopwatch readings. It does not replace
the recorded measurements. It is also not necessarily the truth: it is an
estimate chosen to explain the readings, and the explanation can be wrong.

The process of choosing this row is called **fitting**. We will now work
through the rule ASV uses to decide which fitted row is better.

What a step is
--------------

In our fitted row, the value stays at 10 for a while and then changes to 12.
That change is a **step**. A stretch with one fitted value is a **segment**,
also called a **plateau**. The single value assigned to that stretch is its
**level**.

If we draw a picture with later versions farther to the right and longer
times higher up, our two levels look like this::

    Time in milliseconds

    12 |                  +--------------
       |                  |
    10 | -----------------+
       +---------------------------------
         1   2   3   4    5   6   7   8    Program version

The shape resembles a stair step. **Step detection** means finding where
these lasting changes in level seem to happen.

ASV looks for descriptions made of flat stretches like these. A gradual
change can be harder to describe this way, because it does not have one
clear jump. This is an assumption about the shape of the history, not a fact
that every program's performance must obey.

How to measure how far a fitted value is from a reading
-------------------------------------------------------

Take one reading, 10.1, and a proposed fitted value, 11. Subtract the proposed
value from the reading::

    10.1 - 11 = -0.9

The minus sign says the reading is below the proposal. For now, we only care
about the distance between them. The distance is 0.9 milliseconds.

For another reading, 12.1, the calculation is::

    12.1 - 11 = 1.1

This reading is above the proposal. Its distance from the proposal is 1.1
milliseconds. In both cases, distance is zero or positive.

The mathematical name for taking the size of a number without its minus sign
is **absolute value**. In Python, a programming language, it is written
``abs``::

    abs(-0.9) = 0.9
    abs(1.1)  = 1.1
    abs(0)    = 0

The parentheses tell us which number the instruction applies to. Thus
``abs(10.1 - 11)`` says: first subtract 11 from 10.1, then take the absolute
value of the result. The answer is 0.9.

If you see ``abs(x - y)``, the letters x and y are placeholders for two
numbers. The expression means their distance apart. Reversing the order
does not change that distance: ``abs(10 - 11)`` and ``abs(11 - 10)`` both
equal 1. In mathematical writing you may also see ``|x - y|``; the vertical
bars mean the same thing here.

We will call this distance the **fitting error** for one reading. Error here
does not mean that the program crashed or the measurement was recorded
incorrectly. It just means a difference between the reading and our proposed
description.

Why the one level description has a total error of 8
----------------------------------------------------

Consider a different explanation of our eight measurements: perhaps the
usual timing was 11 milliseconds throughout, with variation above and below
it. For the moment, treat 11 as a proposal to evaluate. A later section
explains where that number comes from.

This explanation has one segment and no change in level::

    Version          1      2      3      4      5      6      7      8
    Measured time   10     10.1    9.9   10     12     12.1   11.9   12
    Proposed time   11     11     11     11     11     11     11     11

Here is every subtraction and every absolute value::

    Version    Measured minus proposed       Distance
       1           10   - 11 = -1.0              1.0
       2           10.1 - 11 = -0.9              0.9
       3            9.9 - 11 = -1.1              1.1
       4           10   - 11 = -1.0              1.0
       5           12   - 11 =  1.0              1.0
       6           12.1 - 11 =  1.1              1.1
       7           11.9 - 11 =  0.9              0.9
       8           12   - 11 =  1.0              1.0

Add the distances for the first four versions::

    1.0 + 0.9 + 1.1 + 1.0 = 4.0

Add the distances for the last four versions::

    1.0 + 1.1 + 0.9 + 1.0 = 4.0

Then add those two totals::

    4.0 + 4.0 = 8.0

This is where 8 comes from. It is the sum of eight separate distances from
the proposed timing. It is not a measurement saying one version took
8 milliseconds. It is not an 8 percent slowdown. It is a total mismatch
used to judge this proposed description.

Why we add distances instead of signed differences
--------------------------------------------------

Imagine just two readings, 10 and 12, with a proposed level of 11. The signed
differences would be -1 and +1. Adding them gives zero::

    -1 + 1 = 0

That zero would hide the fact that both readings differ from the proposal.
Using distances prevents one difference from canceling another::

    abs(-1) + abs(1) = 1 + 1 = 2

That is the reason for ``abs`` in the fitting rule.

Why 11 is a reasonable level for the one segment proposal
---------------------------------------------------------

We chose 11 above to demonstrate the calculation. ASV needs a systematic
way to choose the level for any segment, including ones with thousands of
readings. When each reading counts equally, the rule is to choose a **median**.

To find a median, put the numbers in order from smallest to largest. If there
is an odd number of readings, choose the middle one. For example, the sorted
list ``9, 10, 100`` has middle value 10.

For an even number of readings, the usual convention is to take the two
middle numbers, add them, and divide by two. Our eight readings, sorted, are::

    9.9, 10, 10, 10.1, 11.9, 12, 12, 12.1
                 middle pair: 10.1 and 11.9

There are four numbers on each side of the division between 10.1 and 11.9.
The median calculation is::

    10.1 + 11.9 = 22
    22 divided by 2 = 11

Why use a median? A median gives the smallest possible sum of absolute
distances for a single level. Sometimes more than one level ties for that
smallest total. In this example, any level from 10.1 through 11.9 gives
total error 8; 11 is the middle choice ASV returns.

You may know a different summary called the **mean**, often called the
average: add all the readings and divide by how many there are. A mean and
a median need not be equal. For readings ``10, 10, 100``, the mean is 40,
because ``10 + 10 + 100 = 120`` and 120 divided by 3 is 40. The median is 10.

The total distance from 10 is ``0 + 0 + 90 = 90``. The total distance from
40 is ``30 + 30 + 60 = 120``. The median gives the smaller total distance.
This also shows why an unusually large reading has less influence on a
median than on a mean.

How two segments reduce the total error to 0.4
----------------------------------------------

Return to the explanation that versions 1 through 4 have level 10, and
versions 5 through 8 have level 12. The first group's sorted readings are
``9.9, 10, 10, 10.1``; their middle pair is 10 and 10. Adding them gives 20,
and dividing by two gives median 10. The second group's sorted readings are
``11.9, 12, 12, 12.1``. Its middle pair is 12 and 12: their sum is 24,
and 24 divided by two gives median 12.

Here is every distance for that proposal::

    Version    Measured minus proposed       Distance
       1           10   - 10 =  0.0              0.0
       2           10.1 - 10 =  0.1              0.1
       3            9.9 - 10 = -0.1              0.1
       4           10   - 10 =  0.0              0.0
       5           12   - 12 =  0.0              0.0
       6           12.1 - 12 =  0.1              0.1
       7           11.9 - 12 = -0.1              0.1
       8           12   - 12 =  0.0              0.0

Add the first group's distances::

    0.0 + 0.1 + 0.1 + 0.0 = 0.2

Add the second group's distances::

    0.0 + 0.1 + 0.1 + 0.0 = 0.2

The total is ``0.2 + 0.2 = 0.4``. The two-segment explanation therefore
reduces the mismatch from 8 to 0.4, a reduction of ``8 - 0.4 = 7.6``.

Why the smallest error is not enough
------------------------------------

We could reduce the error even further by copying every reading exactly::

    Version          1      2      3      4      5      6      7      8
    Measured time   10     10.1    9.9   10     12     12.1   11.9   12
    Proposed time   10     10.1    9.9   10     12     12.1   11.9   12

Every reading minus its proposed value is now zero, so the total error is
zero. But the description says that the underlying timing changes at every
version, including each tiny rise and fall of 0.1 milliseconds.

This is called **overfitting**: the explanation follows the particular
measurements so closely that it may be treating ordinary variation as a
meaningful change. If we repeated the experiment, those tiny variations might
look different.

We want a rule that rewards explaining the measurements while discouraging
unnecessary changes. It must balance both aims.

A penalty for adding a change
-----------------------------

ASV adds a numerical **penalty** for each change in fitted level. A penalty
is just an extra amount added to the score; it is not money or extra time
actually spent by the program.

For an example, choose a penalty of 1 per change. Compare three proposals::

    Proposal                   Total error    Changes    Total score
    One level throughout             8           0       8 + 0 = 8
    Level 10, then level 12           0.4         1       0.4 + 1 = 1.4
    Copy every reading               0           7       0 + 7 = 7

There are seven changes in the last proposal because eight versions have
seven gaps between them, and the proposed level changes at every gap.

Smaller scores are preferred. Of these proposals, the two-level description
wins with 1.4. It removes most of the mismatch without claiming a lasting
change for every tiny fluctuation. ASV's exact fitting algorithm also selects
this proposal when it searches all possible divisions for this example and
this penalty.

In the source code, the per-change penalty is called **gamma**. Gamma is just
a parameter name. We will use it now that the meaning of the penalty is clear.

How changing the penalty changes the answer
-------------------------------------------

If gamma is 10 instead of 1, the same proposals have these scores::

    One level:          8 + no changes                    = 8
    Two levels:         0.4 + one change costing 10        = 10.4
    Copy each reading:  0 + seven changes costing 10 each  = 70

Now the one-level description wins. Adding the step saves 7.6 in fitting
error, but costs 10 in the penalty. The saving is not enough to justify it.

At the other extreme, choose a penalty of 0.01. The one-level score stays 8
because there are no changes. The two-level score is ``0.4 + 0.01 = 0.41``.
Copying every reading has zero error and seven changes, so its score is
seven times 0.01, which is 0.07. Changes have become so inexpensive that
following the small fluctuations is cheaper than the two-level description.

This is why choosing gamma matters. An exact search can find the best answer
for a specified gamma, but the best answer under that rule is not automatically
the best explanation of real performance.

Where the Potts model fits into this story
------------------------------------------

The combination we have just built has a mathematical name: a **Potts**
fitting problem. Its defining idea here is to prefer neighboring fitted
values to stay the same, while allowing a change when it improves the fit
enough to pay for the penalty.

The name comes from a model in physics, but ASV does not need to simulate
physics to use this idea. It compares possible descriptions using our score::

    Total score = total fitting error + penalty for all changes

When each reading counts equally, a more compact way to write the same
instruction is::

    Total score = sum of abs(measured time - fitted time)
                  + gamma * number of changes

The word ``sum`` means add up the distances for all readings. The symbol
``*`` means multiply. For the two-level proposal with gamma 1, the second
line contributes ``1 * 1 = 1``.

The charge counts changes, not their sizes. A change from 10 to 12 pays one
gamma, and a change from 10 to 100 also pays one gamma. Their effects on
fitting error differ, so the measurements still affect whether each is useful.

How ASV finds where to put the changes
--------------------------------------

Until now, we proposed the split between versions 4 and 5 ourselves. ASV
must discover a useful split from the measurements. It can consider a split
after version 1, after version 2, and so on. It can also consider multiple
splits or no split.

For every proposed group, it needs a level and the error around that level.
It chooses a median and calculates the distances we have already worked
through. It then compares descriptions using their error and penalties.

Trying every full description separately becomes expensive for long
histories. ASV has an exact method called **dynamic programming** that reuses
the best answers for shorter beginnings of the history. Its basic question
is: if we put the last segment here, what was the best way to describe
everything before it? Reusing those earlier answers avoids starting over
for every proposal.

The normal ASV detection path uses an **approximation** for speed. Approximation
means it may not find the very smallest possible score. It starts with short
segments, joins neighbors when that helps the score, and tries moving their
boundaries. This usually reduces the work, but it can miss a better division.

What the C++ module does
------------------------

Python and C++ are programming languages. Much of ASV is written in Python.
It also has a C++ component for calculations that are repeated many times.

That component helps answer questions such as: what is the median of these
readings, how large is their total error, and which division has the lowest
score for this penalty? It saves answers for groups it has already examined
so it can reuse them. Saving answers this way is called **caching**.

C++ makes those calculations faster. It does not supply knowledge of which
software changes are important. The surrounding Python code still chooses
which penalties to try and which fitted steps to report as regressions.

How ASV chooses the penalty
---------------------------

In ordinary use, you do not have to pick gamma separately for every history.
ASV tries several penalties. Each trial produces a proposed description.
ASV then uses a second scoring rule to choose among those descriptions.

There are two questions, answered in sequence:

1. For this particular penalty, which description fits best?
2. Among the descriptions obtained with different penalties, which seems
   to balance useful detail and ordinary measurement variation best?

The second rule considers how many segments the description needs and what
differences remain between the fitted values and the readings. Those
remaining signed differences are called **residuals**. For a reading of
10.1 and a fitted value of 10, the residual is 0.1.

ASV also allows for neighboring residuals to resemble one another. For
example, the computer may be busy during several successive measurements,
making all of them a little slower. This is called **correlated noise**:
the variations have a relationship instead of acting independently.

The second score uses a mathematical operation called a **logarithm**,
written ``log``. For positive numbers, a larger input gives a larger result,
but the size of a change depends on the ratio of the inputs. For example,
halving a number from 10 to 5
changes its logarithm by the same amount as halving a number from 2 to 1.
This helps the rule compare relative reductions in remaining error. The
implementation reference explains why this operation appears in the score.

A perfect copy of every reading has zero remaining error. Because the
logarithm cannot take zero as a finite input, the score adds a small positive
amount first. That amount is called a **noise floor**. Its purpose is to
keep a perfect fit from receiving an unlimited advantage.

These choices make automatic penalty selection a practical rule, not a
promise that the chosen description is correct. That is why improvements
must be tested on many histories, rather than judged by one attractive fit.

What happens when some measurements are more reliable
-----------------------------------------------------

Our calculations so far treated every reading equally. In practice, repeated
timings for one version might be tightly clustered, while timings for another
version vary widely. ASV uses available information about that uncertainty
to give measurements different influence.

This influence is called a **weight**. A weight multiplies the distance for
a reading before the distances are added. For example::

    Reading A: distance 0.2, weight 1  -> 0.2 multiplied by 1 = 0.2
    Reading B: distance 0.2, weight 3  -> 0.2 multiplied by 3 = 0.6

The same mismatch contributes more to the score for reading B. The fitting
rule therefore has more reason to stay close to B. This makes sense if B
is judged more precise, though that judgment can itself be imperfect.

ASV uses uncertainty estimates from repeated measurements to obtain relative
weights. With different weights, it chooses a **weighted median**, which
finds the middle by accumulated weight instead of just the number of readings.
The central idea remains the same: choose one representative level for each
segment and add up its weighted distances from the readings.

When a fitted step becomes a reported regression
------------------------------------------------

Finding a step and deciding to report a regression are separate decisions.
ASV first estimates the levels. It then examines upward steps and decides
which slowdowns should be reported, considering their size, variation within
the segments, and what happened afterward.

Our example changes from 10 to 12 milliseconds. The absolute increase is::

    12 - 10 = 2 milliseconds

To express that increase as a percentage of the earlier timing, divide the
increase by the earlier timing and multiply by 100::

    2 divided by 10 = 0.2
    0.2 multiplied by 100 = 20 percent

Suppose the reporting threshold is 5 percent. A threshold is a cutoff used
to decide whether a change is large enough to report. Five percent of the
earlier time is::

    5 divided by 100 = 0.05
    10 multiplied by 0.05 = 0.5 milliseconds

The increase of 2 milliseconds is larger than 0.5. In this simple example,
it is also much larger than the variation around either fitted level: each
segment's total error is 0.2 across four readings, or 0.05 per reading. ASV
reports the slowdown between our versions 4 and 5.

Later recovery matters too. Imagine several readings at each of these levels::

    10 milliseconds, then 12, then back to 10

There was a temporary slowdown, but the latest performance has recovered.
Under the reporting rules, this need not remain a regression in the report.
By contrast, going from 10 to 14 and then to 12 still leaves a lasting
slowdown. ASV can report the remaining change from 10 to 12.

The reporting threshold and gamma have different jobs. Gamma helps choose
the fitted description. The reporting threshold helps decide which changes
in that description matter enough to display.

How all the pieces connect
--------------------------

Here is the whole process in order, using the terms we have now introduced::

    Run the same task on different program versions
        |
        v
    Record how long each version takes
        |
        v
    Put the recorded timings in version order
        |
        v
    Try a penalty for adding changes
        |
        v
    Find a description with flat segments
    using distances from the readings plus change penalties
        |
        v
    Repeat with other penalties and choose a description
        |
        v
    Examine its upward changes and later recovery
        |
        v
    Report the slowdowns that pass the reporting rules

The C++ component helps with the repeated fitting calculations in the middle.
It does not run all these stages by itself.

The estimated change is associated with versions of the program. If some
versions have no measurements, ASV may only know that the change happened
somewhere between two measured versions. It cannot recover missing timing
information simply by fitting the available readings.

What we could improve
---------------------

There are three different kinds of improvement to investigate:

* Find a better description for the same penalty, without making the
  calculation too slow. This improves the fitting algorithm.
* Choose penalties and treat measurement variation in a way that finds real
  changes more reliably. This improves the selection rule and noise assumptions.
* Make sure the report highlights useful lasting slowdowns. This improves
  the final reporting rules.

To test an idea, we can make up a history where we know when the usual timing
changes, add small random variations, and see whether the detector finds the
change. We must also test histories where the usual timing never changes:
a method that always declares a slowdown would do badly on those histories.

A lower mathematical score is not proof of a more truthful explanation.
It only says that the description is better under that particular scoring
rule. Good experiments check both the calculation and the usefulness of
the resulting reports.

Where to go next
----------------

The `improvement plan <improvement_plan.rst>`_ starts with a plain-language
description of the proposed experiments, then gives the detailed research plan.

The `implementation reference <implementation_details.rst>`_ connects the
ideas here to ASV function names, the exact mathematical formulas, and the
C++ module. It includes executable examples. Its list positions start at zero
because Python numbers list entries that way: version 1 in this guide is
position 0 there, and the change between versions 4 and 5 is between positions
3 and 4. The measurements and calculations are otherwise the same.

These explanations describe source revision
``d33754e129c025beb5c2ca440c3c280433b264f7``. The reference links to the relevant
source files and the original ASV documentation for readers ready to explore
the implementation.
