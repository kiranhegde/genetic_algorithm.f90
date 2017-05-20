!---------------------------------------------------
! PHY2071 – Genetic Algorithms
! URN 6317482
! Name: Anthony Hills
! Jan 26th 2017
!---------------------------------------------------

! This program applies a genetic algorithm to find the optimal solution for
! thrust and angle for a projectile to reach a user-specified target, located
! in the x axis.

PROGRAM genetic_algorithm
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: max_pop=1000, n_max=1000
  REAL, PARAMETER :: pi=2*asin(1.0) !,  tx_min=100.0, tx_max=500.0
  REAL, PARAMETER :: angle_min=0, angle_max=(pi/2), thrust_min=10.0, thrust_max=250.0
  REAL, PARAMETER :: mutation_chance=0.33, crossover_chance=0.33
  REAL, PARAMETER :: accuracy = 0.001
  INTEGER :: n, c, tournament_size, done, t
  REAL :: xs, tx
  REAL, DIMENSION(1:2) :: parent1, parent2, child1, child2, mutated_child
  REAL, DIMENSION(1:2) :: fittest_individual, weakest_individual
  REAL, DIMENSION(1:max_pop, 1:2) :: current_generation, next_gen, random_population
  REAL :: crossover_val, mutation_val
  INTEGER, PARAMETER :: sample_size=10
  INTEGER :: s, e
  REAL, DIMENSION(1:sample_size) :: sample
  REAL :: sum_iterations, mean_generations, the_sum, std_generations
  REAL :: p = 0.0
  INTEGER, PARAMETER :: e_replicants = 1


  WRITE(6,*) "Please enter the x coordinates (value between 100 and 500)"
  WRITE(6,*) "for the target location:"
  READ(5,*) tx
  WRITE(6,*) ""

!==============================================================================


!DO t=1,100
!WRITE(6,*) t
!tournament_size = t
!sum_iterations=0.0
!mean_generations=0.0
!the_sum=0.0
!std_generations=0.0


OPEN(unit=30,file='stats.dat')



DO s=1, sample_size
  WRITE(6,*) s,"/",sample_size

OPEN(unit=40, file='fitness_against_generation_number.dat')

  done=0
  n = 1

  CALL generate_random_population(max_pop, angle_min, angle_max, &
  & thrust_min, thrust_max, random_population)
  
  current_generation = random_population
  ! Sets the current_generation as the randomly generated population

  fittest_individual = current_generation(1,:)
  ! Arbitrary initialization

  DO WHILE ((fitness(fittest_individual,tx) >= accuracy) .AND. (n<=n_max))
  ! Loops until either the fittest individual is evaluated within specified accuracy,
  ! or until the maximum number of generation iterations has been reached.

OPEN(unit=20,file='fitness.dat')
    WRITE(20,*) s,n,fitness(fittest_individual,tx)

    c = 1
    ! Denotes the array position for the individuals within the next_gen array

    CALL evaluate_fittest(current_generation, max_pop, fittest_individual, tx)
    
!DO e=1, e_replicants
next_gen(c,:) = fittest_individual
    ! This is elitism.
    ! Adds the fittest individual from the current generation to the next_gen array.
    c=c+1
!END DO

    DO WHILE (c<=max_pop)
    ! Fills the rest of the next_gen array using selected individuals
    ! Loops until the next_gen array is full

      tournament_size = 4

      CALL tournament(max_pop, tournament_size, current_generation, &
      & parent1, parent2, tx)
      CALL RANDOM_NUMBER(crossover_val)
      CALL RANDOM_NUMBER(mutation_val)

      IF (crossover_val <= crossover_chance) THEN
        CALL crossover(parent1, parent2, child1, child2)
        next_gen(c,:) = child1
        c = c + 1
        next_gen(c,:) = child2
        c = c + 1
      END IF
      IF (mutation_val <= mutation_chance) THEN
        CALL mutate(parent1, mutated_child,n)
        next_gen(c,:) = mutated_child
        c = c + 1
      END IF
      
      IF ((crossover_val > crossover_chance) &
      &  .AND. (mutation_val > mutation_chance)) THEN
        next_gen(c,:) = parent1
        c = c + 1
        next_gen(c,:) = parent2
        c = c + 1
      END IF
    END DO
    ! The next_gen array is now full

    current_generation = next_gen
    ! Sets the current_generation to the newly evolved generation
    n=n+1
    ! Updates the generation number
    CALL evaluate_fittest(current_generation, max_pop, fittest_individual, tx)
    ! Uses the most recent generation to find the newest fittest individual

      WRITE(40,*) n, mean_fitness(current_generation, max_pop, tx), std_fitness(current_generation, max_pop, tx) !, std_fitness

  END DO

	done=1

  CLOSE(40)
  
  sample(s) = n-1
  ! Saves the number of evolution iterations to reach the solution. 

  WRITE(30,*) sample(s) 


END DO

  CLOSE(30)
  CLOSE(20)

  DO s=1, sample_size
    sum_iterations = sum_iterations + sample(s)
  END DO
  mean_generations = sum_iterations/sample_size
  ! Evaluates the mean generations to convergence.

  DO s=1, sample_size
    the_sum = the_sum + (sample(s) - mean_generations) ** 2
  END DO
  std_generations = SQRT(the_sum/(sample_size - 1))
  ! Evaluates the std in the # of generations to convergence.

  WRITE(6,*) "generations to convergence: "
  WRITE(6,*) mean_generations, "+-", std_generations

  CALL evaluate_fittest(current_generation, max_pop, fittest_individual, tx)




! CLOSE(40)

!==============================================================================


!==============================================================================
CONTAINS
  SUBROUTINE generate_random_population(max_pop, angle_min, angle_max, &
    & thrust_min, thrust_max, random_population)
  ! Generates a random population.
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: max_pop
    REAL, INTENT(IN) :: angle_min, angle_max, thrust_min, thrust_max
    REAL, DIMENSION(1:max_pop, 1:2), INTENT(OUT) :: random_population
    REAL, DIMENSION(1:2) :: individual
    REAL :: angle, thrust, r
    INTEGER :: i

    DO i=1,max_pop
      CALL RANDOM_NUMBER(r)
      angle = (r * (angle_max-angle_min)) + angle_min
      CALL RANDOM_NUMBER(r)
      thrust = (r * (thrust_max-thrust_min)) + thrust_min
      ! Generates random values for angle and thrust between their stated
      ! maxmum and minimum value conditions.
      individual(1) = angle
      individual(2) = thrust
      ! Assigns these generated angle and thrust values to the individual.
      random_population(i,:) = individual
      ! Adds the randomly generated individual to the random population. 
    END DO
  END SUBROUTINE


  REAL FUNCTION fitness(individual, tx)
  ! Takes in the value of the individual, and returns its fitness.
  ! This value is a representation of how close this solution is to the
  ! best solution.
    IMPLICIT NONE
    REAL, DIMENSION(1:2), INTENT(IN) :: individual
    REAL, INTENT(IN) :: tx

    CALL launch_projectile(individual, xs, done)
    fitness = ABS(tx-xs)
    ! Low fitness approaches the solution.
  END FUNCTION fitness


  REAL FUNCTION mean_fitness(population, population_size, tx)
  ! Returns the mean fitness of the given population.    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: population_size
    REAL, INTENT(IN) :: tx
    REAL, DIMENSION(1:population_size, 1:2), INTENT(IN) :: population
    REAL :: sum_fitness
    INTEGER :: i

    sum_fitness = 0

    DO i=1, population_size
      sum_fitness = sum_fitness + fitness(population(i,:), tx)
    END DO

    mean_fitness = sum_fitness/ i
  END FUNCTION mean_fitness


  REAL FUNCTION std_fitness(population, population_size, tx)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: population_size
    REAL, DIMENSION(1:population_size, 1:2), INTENT(IN) :: population
    REAL, INTENT(IN) :: tx
    REAL :: mean
    REAL :: the_sum
    INTEGER :: i    

    mean = mean_fitness(population, population_size, tx)

    DO i=1, population_size
      the_sum = the_sum + (fitness(population(i,:), tx) - mean) ** 2   
    END DO

    std_fitness = SQRT(the_sum/ (population_size - 1))
  END FUNCTION std_fitness


  SUBROUTINE evaluate_fittest(some_population, population_size, fittest_individual, tx)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: population_size
    REAL, DIMENSION(1:population_size, 1:2), INTENT(IN) :: some_population
    REAL, INTENT(IN) :: tx
    REAL, DIMENSION(1:2), INTENT(OUT) :: fittest_individual
    REAL :: last_fitness
    INTEGER :: k

    !last_fitness=10000.00
    ! Arbitrary initialization
    fittest_individual = some_population(k,:)
    ! Arbitrary initialization 

    DO k=1,population_size
    ! This loop finds the fittest individual from the population
      IF (fitness(some_population(k,:),tx) <= & 
         & fitness(fittest_individual,tx)) THEN
        fittest_individual = some_population(k,:)
      END IF
      ! last_fitness = fitness(some_population(k,:),tx)
    END DO
  END SUBROUTINE


  SUBROUTINE evaluate_weakest(some_population, population_size, weakest_individual, tx)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: population_size
    REAL, DIMENSION(1:population_size, 1:2), INTENT(IN) :: some_population
    REAL, INTENT(IN) :: tx
    REAL, DIMENSION(1:2), INTENT(OUT) :: weakest_individual
    REAL :: last_fitness
    INTEGER :: k

    last_fitness=-1.0
    ! Arbitrary initialization

    DO k=1,population_size
    ! This loop finds the weakest (i.e. largest value in fitness) 
    ! individual from the population
      IF (fitness(some_population(k,:),tx) >= last_fitness) THEN
        weakest_individual = some_population(k,:)
      END IF
      last_fitness = fitness(some_population(k,:),tx)
    END DO
  END SUBROUTINE

  INTEGER FUNCTION array_position_of(individual, population, population_size)
  ! Returns the position of an individual in the array of a given population.
    IMPLICIT NONE  
    INTEGER, INTENT(IN) :: population_size
    REAL, DIMENSION(1:2), INTENT(IN) :: individual
    REAL, DIMENSION(1:population_size, 1:2), INTENT(IN) :: population
    INTEGER :: i

    DO i=1, population_size
      IF ((population(i,1) == individual(1)) .AND. (population(i,2) == individual(2))) THEN
      ! Checks if the individuals are the same.
        array_position_of = i
        ! Saves the position of the found individual in the array. 
        EXIT
      END IF
    END DO
  END FUNCTION


  SUBROUTINE elitism(max_pop, old_population, new_population, tx)
  ! Ensures that the fittest individual from a population is added to the
  ! next generation, i.e. the new population.    
   IMPLICIT NONE
    INTEGER, INTENT(IN) :: max_pop
    REAL, INTENT(IN) :: tx
    REAL, DIMENSION(1:max_pop, 1:2), INTENT(IN) :: old_population    
    REAL, DIMENSION(1:max_pop, 1:2), INTENT(OUT) :: new_population
    REAL, DIMENSION(1:2) :: fittest_individual, weakest_individual 
    INTEGER :: a, b

    CALL evaluate_fittest(old_population, max_pop, fittest_individual, tx)
    CALL evaluate_weakest(new_population, max_pop, weakest_individual, tx)

    a = array_position_of(weakest_individual, new_population, max_pop)
    b = array_position_of(fittest_individual, old_population, max_pop)

    new_population(a,:) = old_population(b,:)   

    ! Replaces the weakest_individual in the new population with the fittest of
    ! the old population. 
  END SUBROUTINE


  SUBROUTINE tournament(max_pop, tournament_size, current_generation, &
    & parent1, parent2, tx)
  ! Generates two successful parents out of tournament selection.
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: tournament_size, max_pop
    REAL, INTENT(IN) :: tx
    REAL, DIMENSION(1:max_pop, 1:2), INTENT(IN) :: current_generation
    REAL, DIMENSION(1:2), INTENT(OUT) :: parent1, parent2
    REAL, DIMENSION(1:tournament_size, 1:2) :: tournament_population
    REAL, DIMENSION(1:2, 1:2) :: fittest_individuals
    INTEGER :: p, i, random_integer
    REAL :: r
    ! If the tournament size is larger, weak individuals have a smaller chance
    ! to be selected.

    DO p=1,2
    ! This loop finds two parents for reproduction.
      DO i=1, tournament_size
      ! This loop fills the tournament population array.
        CALL RANDOM_NUMBER(r)
        random_integer = CEILING(SIZE(current_generation)*r)
        ! Creates a random integer from 1 to the population size.     
        tournament_population(i,:) = current_generation(random_integer,:)
        ! Picks at random, an individual from the current_generation and adds
        ! it to the tournament.
      END DO    

      CALL evaluate_fittest(tournament_population, tournament_size, fittest_individual, tx)
      fittest_individuals(p,:) = fittest_individual
    END DO

    parent1 = fittest_individuals(1,:)
    parent2 = fittest_individuals(2,:)
    ! Acquires two successful individuals from two seperate tournaments
    ! and prepares them for reproduction.
  END SUBROUTINE tournament


  SUBROUTINE crossover(parent1, parent2, child1, child2)
  ! Two parents produce offspring via crossover.
    REAL, DIMENSION(1:2), INTENT(IN) :: parent1, parent2
    REAL, DIMENSION(1:2), INTENT(OUT) :: child1, child2

    child1(1) = parent2(1)
    child1(2) = parent1(2)
    ! Crosses the parents first set of possible values for angle and 
    ! thrust to produce the first child.

    child2(1) = parent1(1)
    child2(2) = parent2(2)
    ! Crosses the parents second possible set of values for angle and
    ! thrust to produce a second child. 
  END SUBROUTINE crossover


  SUBROUTINE mutate(parent, mutated_child, n)
    REAL, DIMENSION(1:2), INTENT(IN) :: parent
	INTEGER, INTENT(IN) :: n
    REAL, DIMENSION(1:2), INTENT(OUT) :: mutated_child
    REAL :: delta_angle, delta_thrust,delta_tamp,delta_tangle
	REAL :: r
    ! Perhaps have delta_angle and delta_thrust as intent in variables.
    ! Adjust these delta values as functions, decreasing in terms of percentage
    ! relative to angle and thrust size.

	CALL RANDOM_NUMBER(r)
	delta_tangle=EXP(-(REAL(n))/15)
    delta_angle = delta_tangle*(r-0.5)
	CALL RANDOM_NUMBER(r)
    delta_tamp=EXP(-(REAL(n))/200)
    delta_thrust = delta_tamp*(r-0.5)

    mutated_child(1) = parent(1) + delta_angle
    mutated_child(2) = parent(2) + delta_thrust
    ! Adds a small random value, delta, within an accuracy to the child.
    ! Ensure that the mutated child value does not exceed BCs.

    ! Perhaps have delta large for low n, and as n increases, narrow delta.
    ! Do this in a logarithmic fashion.
  
    !WRITE(6,*) "angle, then thrust"
   ! WRITE(6,*) parent(1), delta_angle
  !  WRITE(6,*) parent(2), delta_thrust
 !   WRITE(6,*) "end mutation"
  
  END SUBROUTINE mutate

  REAL FUNCTION reduce_angle(angle)
  ! Returns the smallest angle in radians, instead of the angle 
  ! calculated through multiple revolutions. 
    IMPLICIT NONE
    REAL, PARAMETER :: pi=3.141592653589793238462643383279502884197169
    REAL, INTENT(IN) :: angle ! Units in rad.
    INTEGER :: n

    n = 0
    ! n corresponds to the number of revolutions through the circle.

    DO WHILE (69 == 69)
    ! i.e. forever
      n = n + 1
      IF (2*pi*n - ABS(angle) >= 0) THEN  
        EXIT
        ! Exits once there are no more revolutions possible.
      END IF
    END DO

    n = n - 1
    reduce_angle = angle - 2*pi*n
  END FUNCTION reduce_angle

!______________________________________________________________________________

  SUBROUTINE launch_projectile(individual, xs, done)
  ! Evaluates the location on the x-axis that the projectile lands, xs,
  ! when given values for the angle and thrust.
    IMPLICIT NONE
    REAL, DIMENSION(1:2), INTENT(IN) :: individual
	  INTEGER, INTENT(IN) :: done
    ! Input for the solution
    REAL, INTENT(OUT) :: xs
    ! Output : the x position where the projectile hits the ground (y=0)
    REAL :: x, y, vx, vy, ax, ay
    ! Variable used for the integration. Self explanatory
    REAL :: lx, ly
    ! Previous position in the integration. We keep track of this
    ! to interpolate the last timestep. The integration will more
    ! than certainly give us a last x where y < 0. We interpolate
    ! to find the x position where y=0
    REAL, PARAMETER :: g        = -9.81
    REAL, PARAMETER :: friction = 0.05
    REAL, PARAMETER :: dt       = 0.01
    ! Constants.
    REAL :: thrust, angle

    angle = individual(1)
    thrust = individual(2)


    ! Initialization :
    !  - last and current positions set to (0, 0)
    !  - velocity set to the initial conditions
    lx = 0.0
    ly = 0.0
    x = 0.0
    y = 0.0
    vx = COS(angle) * thrust
    vy = SIN(angle) * thrust

	IF (done==1) THEN
		  OPEN(unit=20,file='trajectory.dat')
      WRITE(20,*) x,y
	END IF

    ! Looping until hitting the ground
    DO
     ! Updating last position
     lx = x
     ly = y

     ! Computing the acceleration. On x, only friction, on y gravity.
     ax = -vx * 0.05
     ay = g

     ! Updating velocities, first order Euler
     vx = vx + ax * dt
     vy = vy + ay * dt

     ! Updating positions, second order Euler
     x = x + vx * dt + ax * 0.5 * dt**2
     y = y + vy * dt + ay * 0.5 * dt**2

	IF (done==1) THEN
		  OPEN(unit=20,file='trajectory.dat')
      WRITE(20,*) x,y
	END IF

     ! If we hit the ground, we exit the loop
     IF (y <= 0.0) EXIT
    END DO

		IF (done==1) THEN
	CLOSE(20)
	END IF

    ! Interpolating the result to get the intersection with the origin
    xs = lx + (x - lx) * ly / (ly - y)
  END SUBROUTINE launch_projectile
END PROGRAM genetic_algorithm
