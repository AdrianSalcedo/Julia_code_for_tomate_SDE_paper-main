function Verify_R0d_greatherthan_one_condition(par)
   #Input:
   #   par: vector with parameters
   # Output:
   #   deterministic R0: deterministic squart reproductive number

   Rd0 = Compute_deterministic_R0(par)

   cond = 1.5 > Rd0 > 1.0
   #test = sign(cond) == 1.0
   return cond #test
end
