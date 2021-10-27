open Complex

let nextpowof2 n =
   let x = n in 
   let x = x lor (Int.shift_right x 1) in 
   let x = x lor (Int.shift_right x 2) in
   let x = x lor (Int.shift_right x 4) in 
   let x = x lor (Int.shift_right x 8) in 
   let x = x lor (Int.shift_right x 16) in 
   x + 1
  
  
let rec fft a n w =
   if w = one || n = 1 then a else
      let ae = Array.init (n / 2) (fun i -> a.(2 * i)) in
      let ao = Array.init (n / 2) (fun i -> a.(2 * i + 1)) in
      let r1 = fft ae (n / 2) (mul w w) in
      let r2 = fft ao (n / 2) (mul w w) in
      for i = 0 to (n / 2) - 1 do
         a.(i) <- add r1.(i) ( mul ( pow w {re = float_of_int i; im = 0.0} ) r2.(i) );
         a.(i + (n / 2)) <- sub r1.(i) ( mul ( pow w {re = float_of_int i; im = 0.0} ) r2.(i) )
      done;
      a

let multiply_pol () = 

   (* Reading the inputs *)
   print_endline "Degree of first polynomial :";
   let n1 = read_int() in 
   print_endline "Degree of second polynomial :";
   let n2 = read_int() in 

   (* Evaluating the degree of the product polynomial *)
   let n = nextpowof2 (n1 + n2) in

   let a1 = Array.make n { re = 0.0; im = 0.0 } in
   let b1 = Array.make n { re = 0.0; im = 0.0 } in

   (* First polynomial *)
   print_endline "Coefficient vector of 1st polynomial :";
   for i = 0 to n1 do 
      a1.(i) <- {re = read_float(); im = 0.0 };
   done;

   (* Second polynomial *)
   print_endline "Coefficient vector of 2nd polynomial :";
   for i = 0 to n2 do 
      b1.(i) <- {re = read_float(); im = 0.0 };
   done;
    
   (* Evaluating omega w *)
   let theta = 2.0 *. Float.pi /. float_of_int n in
   let w = {re = cos theta; im = sin theta } in 

   (* FFT of polynomial gives Value Repn of each polynomial *)
   let a2 = fft a1 n w in
   let b2 = fft b1 n w in 
   
   (* Evaluating Product polynomial in Value Repn *)
   let c1 = Array.init n (fun i -> mul a2.(i) b2.(i)) in

   (* Interpolation : Value -> Coeff *)
   let win = conj w in 
   let c2 = fft c1 n win in
   let d = Array.map (fun c -> { re = c.re /. (float_of_int n); im = 0. }) c2 in

   (* Printing the Product *)
   print_endline "\nCoefficient vector of product :" ;
   Array.iteri (fun i di -> if i < (n1 + n2 + 1) then Printf.printf "%.0f " di.re) d;
   print_newline ();
   print_newline ();

   print_endline "The product polynomial :";
   Array.iteri 
   (fun i di -> 
      if i < (n1 + n2) then Printf.printf "%.0fx^%d + " di.re i  
      else if i = (n1 + n2) then Printf.printf "%.0fx^%d" di.re i) d;
   print_newline ();
   print_newline ()
   
let () = multiply_pol ()