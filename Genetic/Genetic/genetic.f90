!***************************************************************************
!MAIN GENETIC PROGRAM
!Description: Genetic algorithm to match a randomly generated sequence 
!to an input string. 
!Methods:
!   Scorer    -> Scores string for algorithm to find "Most Fit Candidates"
!   Minimum   -> Finds minimum score
!   CrossOver -> Breeds most fit randomly with population after discarding
!                the least fit members.     
!   Mutation  -> Randomly mutates values during crossOver events
!***************************************************************************
program genetic
    implicit none
    
    !Declare loop variables
    integer i, j 
    
    character wordList(11,1000)
    CHARACTER*11 key
    integer score(1000), min, loop
    integer NumberOfChromosomes, wordLength, value, flag
    real rand
    integer, parameter :: seed = 777
    call srand(seed)
    
    !Initialize variables
    NumberOfChromosomes = 1000
    wordLength = 11
    key = "NASA"
    
    !Generate Random set of character data
    do 10 j = 1, NumberOfChromosomes
        do 20 i = 1, wordLength
            if (rand(0).LT.0.6) then  
                wordList(i,j) = char(int((rand(0)*(122+1-97))+97))
            else if (rand(0).GT.0.6.AND.rand(0).LT.0.95) then  
                wordList(i,j) = char(int((rand(0)*(90+1-65))+65))
            else 
                wordList(i,j) = char(32)
            endif
20      continue
10  continue    
   
    !score the data set
    call scorer(wordList, score)
    
    !initialize min
    call minimum(score, min)

    !Loop over until score of one value of scorer is exactly zero, 
    !which means the strings are the same
    do while(min.GT.0)  
        
        !crossOver the most fit values
        call crossOver(wordList, score)
        
        !mutate some randomly
        call mutation(wordList)
    
        !score the data set
        call scorer(wordList, score)
        
        !Get Best Score
        call minimum(score, min)
        
        !Print status every 100 generations
        if (mod(loop,100).EQ.0) then 
            call system('CLS')
            print 1, "********************"
            print 2, "Generation: ", loop
            print 3, "Minimum Score: ", min
            do j = 1,1
                print *, (wordList(i,j), i = 1,wordLength) 
            enddo
            print 1, "********************"
1           format(A20)            
2           format(A12, I8)
3           format(A15, I5)   
        endif

        !keep track of loops
        loop = loop + 1
    enddo
    
    !Output information 
    print 5, " It took ", loop, " generations to finish"
5   format(A9, I7, A22)  
    do j = 1,1000
        print *, (wordList(i,j), i = 1,wordLength)         
    enddo
    
    end program genetic
    
    !*******************************************************************
    !MINIMUM SUBROUTINE
    !*******************************************************************
    subroutine minimum(score, min)
    integer j, min, score(1000)
    min = 10000
    do j = 1, 1000
        if (min.GT.score(j)) then 
            min = score(j)
        endif
    enddo
    return 
    end subroutine minimum
   
    !*******************************************************************
    !SCORER SUBROUTINE
    !*******************************************************************
    subroutine scorer(wordList, score)
    integer i, j, m, n
    character wordList(11,1000)
    CHARACTER*11 key
    integer score(1000)
    integer NumberOfChromosomes, wordLength, value
    
    !Initialize variables
    NumberOfChromosomes = 1000
    wordLength = 11
    key = "NASA"
    
    !Clear score table between runs
    do j=1, NumberOfChromosomes
        score(j) = 0
    enddo
    
    !score the values 
    do 12 j = 1, NumberOfChromosomes
        do 22 i = 1, wordLength
            !Difference of ascii(key[i]) and ascii(word[i])
            if ((ichar(key(i:i)) - ichar(wordList(i,j))).LT.0) then
                value = (-1)*(ichar(key(i:i)) - ichar(wordList(i,j)))
            else
                value = (ichar(key(i:i)) - ichar(wordList(i,j)))
            endif
            score(j) = score(j) + value
22      continue
12  continue   
    return 
    end subroutine scorer
    
    !********************************************************************
    !CROSSOVER SUBROUTINE
    !Methods:
    !   BoolLeast    -> returns a true/false flag if value is in list
    !   Bool         -> returns a true/false flag if value is in list
    !   Probability  -> calculates the probability of the fit members to 
    !                   breed
    !   CrossingPair -> Determines which fit member is to be crossed with
    !                   a given population member
    !********************************************************************
    subroutine crossOver(wordList, score)

    integer i, j, k, max, index, value, min
    character wordList(11,1000)
    character*11 key
    character mostFit(11,5), leastFit(11,1000)
    integer score(1000), mostFitScore(5), indexList(5), leastFitScore(50)
    integer leastFitIndexList(50)
    integer NumberOfChromosomes, wordLength
    logical flag, bool, boolLeast
    real probabilityList(5), prob, rand, roll
    
    !Initialize variables 
    key = "NASA"
    NumberOfChromosomes = 1000
    wordLength = 11
    min = 10000
    max = 0
    
    !find the best scores and return as a list
    do k = 1, 5
        do j = 1, NumberOfChromosomes
            flag = bool(j, indexList)
            if ((min.GT.score(j)).AND.(flag.NE..TRUE.)) then 
                min = score(j)
                index = j
            endif
        enddo
        mostFitScore(k) = min
        indexList(k) = index
        min = 10000
    enddo
    
    !find the worst scores and return as a list
    do k = 1, 50
        do j = 1, NumberOfChromosomes
            flag = boolLeast(j, leastFitIndexList)
            if ((max.LT.score(j)).AND.(flag.NE..TRUE.)) then 
                max = score(j)
                index = j
            endif
        enddo
        leastFitScore(k) = max
        leastFitIndexList(k) = index
        max = 0
    enddo
    
    !Kill off the least fit and allow new chromosmes to enter the gene pool
    do j = 1, 50
        do i = 1, wordLength
            if (rand(0).LT.0.6) then  
                wordList(i,leastFitIndexList(j)) = char(int((rand(0)*(122+1-97))+97))
            else if (rand(0).GT.0.6.AND.rand(0).LT.0.95) then  
                wordList(i,leastFitIndexList(j)) = char(int((rand(0)*(90+1-65))+65))
            else 
                wordList(i,leastFitIndexList(j)) = char(32)
            endif
        enddo
    enddo
       
    !Build Most Fit List
    do k = 1,5
        do i = 1, wordLength
            index = indexList(k)
            mostFit(i,k) = wordList(i,index)
        enddo
    enddo
    
    !Generate Probability List for crossover events 
    call probability(mostFit, probabilityList, mostFitScore) 
    
    !Cross the substrings to generate a new population
    do j = 1, NumberOfChromosomes
        call crossingPair(probabilityList, value)
        do i = 1, wordLength
            roll = (int(rand(0)*11)+1)
            if (i.LT.roll) then 
                wordList(i,j) = mostFit(i, value)
            endif
            !fix characters that are correct in fit population
            if ((ichar(key(i:i)) - ichar(mostFit(i,value))).EQ.0) then 
                wordList(i,j) = mostFit(i,value)
            endif
        enddo
    enddo
    return 
    end subroutine crossOver
    
    !*******************************************************************
    !BOOLLEAST FUNCTION
    !*******************************************************************
    logical function boolLeast(index, leastFitIndexList) 
    
    !Declare variables
    integer leastFitIndexList(50), index, i 
    
    !Initialize variables
    boolLeast = .FALSE.
    
    do i = 1,50
        if (leastFitIndexList(i).EQ.index.AND.index.NE.0) then
            boolLeast = .TRUE.
        endif
    enddo
    end function boolLeast
    
    
    !*******************************************************************
    !BOOL FUNCTION
    !*******************************************************************
    logical function bool(index, indexList) 
    
    !Declare variables
    integer indexList(5), index, i 
    
    !Initialize variables
    bool = .FALSE.
    
    do i = 1,5
        if (indexList(i).EQ.index.AND.index.NE.0) then
            bool = .TRUE.
        endif
    enddo
    end function bool
    
   
    !*******************************************************************
    !PROBABILITY SUBROUTINE
    !*******************************************************************
    subroutine probability(mostFit, probabilityList, mostFitScore)
    
    !Declare variables
    real sum, probabilityList(5)
    integer i, j, mostFitScore(5)
    character*11 key
    character mostFit(11,5)
    
    !Initialize variables
    sum = 0
    key = "nasa"
    
    !Sum scores for probability
    do i = 1,5
        sum = sum + mostFitScore(i)
    enddo    
    
    !Find the probabilty of the event not occuring and use it as the probability of interest 
    !as the best values are the smallest. 
    do i = 1,5
        probabilityList(i) = 1 - mostFitScore(i)/sum
    enddo
    return 
    end subroutine probability
    
    !*******************************************************************
    !CROSSINGPAIR SUBROUTINE
    !*******************************************************************
    subroutine crossingPair(probabilityList, value)
    
    !Declare variables
    integer k, value, arg1, arg2, range
    real roll, probabilityList(5), rand
    
    !Initialize variables
    roll = (rand(0)*100) + 1
    arg1 = 100
    arg2 = probabilityList(1)
    range = 0
    
    !get value from probabilty table as 0 < roll < range(arg1-arg2).
    !Then update the start and end to range < roll < range(arg2 -arg3) etc.
    !Finally, note the args reference k+1, so if all the previous checks fail, 
    !the defualt must be the end of the array, i.e. 5. and we do not loop over 
    !that value so the range is len(list)-1. 
    do k = 1,4
        if (roll.GT.range.AND.roll.LT.(arg1 - arg2)) then 
            value = k
        else
            value = k+1
        endif
        arg1 = probabilityList(k)
        arg2 = probabilityList(k+1)
    enddo
    
    return
    end subroutine crossingPair
    
    !*******************************************************************
    !MUTATION SUBROUTINE
    !*******************************************************************
    subroutine mutation(wordList)
    character wordList(11,1000)
    real rand, roll
    integer i,j, NumberOfChromosomes
    
    !Initialize the values
    NumberOfChromosomes = 1000
    
    !roll the rand
    roll = rand(0)
    
    !Mutate random genes
    do j = 1, NumberOfChromosomes
        do i = 1, 11
            !Mutation rate
            if (roll < .05) then 
                if (rand(0).LT.0.6) then  
                    wordList(i,j) = char(int((rand(0)*(122+1-97))+97))
                else if (rand(0).GT.0.6.AND.rand(0).LT.0.95) then  
                    wordList(i,j) = char(int((rand(0)*(90+1-65))+65))
                else 
                    wordList(i,j) = char(32)
                endif
            endif
        enddo
    enddo
    return 
    end subroutine mutation
    
    
    
    
   