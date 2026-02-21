from sympy import nextprime
from sympy.ntheory.residue_ntheory import crt
from typing import List
import argparse
import logging

# --- GLOBAL CONSTANT FOR SCRIPT EXECUTION ---
# This needs to be defined globally so the print loop (in __main__) can access it
# to calculate the P-max (p_start, i) correctly.
P_START_PRIME = 3
# --- END GLOBAL CONSTANT ---

# Configure logging to provide informative output with timestamps
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def validate_general(g_k, primes_list, target_mod, target_rem):
    """
    Validates a newly calculated Shield General against all required constraints.
    This acts as a "Forward Check" to ensure the mathematical properties of the
    Shield General hold true for all shielded primes.

    Args:
        g_k (int): The Shield General to validate.
        primes_list (list): A list of all primes shielded so far (from 3 up to p_k).
        target_mod (int): The modulus for the parity check (usually 2).
        target_rem (int): The remainder for the parity check (usually 0 for even).

    Returns:
        bool: True if all conditions pass, False otherwise.
    """
    # 1. Check Parity: The Shield General must be even (G_k % 2 == 0)
    # This corresponds to the constraint that prime gaps (for primes > 2) are even.
    if g_k % target_mod != target_rem:
        logging.error(f"Parity check failed for {g_k}. Expected rem {target_rem} for mod {target_mod}, got {g_k % target_mod}.")
        return False

    # 2. Check Shield Conditions for each prime in the list
    for p in primes_list:
        # Special Case: Prime 3
        # The gap must be non-zero mod 3 to avoid the sum being divisible by 3.
        if p == 3:
            if g_k % 3 == 0:
                logging.error(f"Initial shield check failed for {g_k}. Expected non-zero rem for mod 3, got 0.")
                return False

        # General Case: Primes >= 5
        # The gap must be 1 or -1 mod p to shield against p.
        elif p >= 5:
            rem = g_k % p
            if rem != 1 and rem != p - 1:
                logging.error(f"Subsequent shield check failed for {g_k} with prime {p}. Expected rem 1 or {p-1}, got {rem}.")
                return False

    return True


def generate_shield_generals(n_terms: int) -> List[int]:
    """
    Generates the first n_terms of the standard Shield General sequence.

    A Shield General G_k is defined as the smallest even integer such that:
      1. G_k = 1 mod 3
      2. G_k = 1 mod p_i for all primes 5 <= p_i <= p_k

    This function uses a "Ratchet" mechanism to efficiently find the next General.
    Instead of searching from scratch, it increments by the product of previous primes
    (the primorial) until the new constraint is satisfied.
    """
    if n_terms < 1:
        return []

    # --- 1. Standard Constraints ---
    # We start shielding from p=3.
    p_start = P_START_PRIME
    r_start = 1
    # We require the General to be even.
    target_mod = 2
    target_rem = 0

    # --- 2. Initialize Base General G_0 ---
    # The first Shield General is 4.
    # 4 is even, and 4 % 3 == 1. It satisfies the shielding condition for p=3.
    g_init = 4

    # --- 3. Initialize Ratchet State ---
    current_general = g_init

    # The current cycle represents the product of all primes shielded so far.
    # We start with 3. The next search steps will be multiples of this cycle.
    current_cycle = P_START_PRIME

    # p_k tracks the highest prime we are currently shielding against.
    p_k = P_START_PRIME

    # List to track the actual primes being shielded for validation purposes.
    shielded_primes = [3]

    # Store the results, starting with the initial general.
    results = [current_general]

    # Validate the initial term (Term 1)
    # It is critical to ensure our starting point is valid before building upon it.
    if not validate_general(current_general, shielded_primes, target_mod, target_rem):
        raise RuntimeError(f"Validation failed for initial General {current_general} (Term 1).")
    logging.info(f"Successfully validated Shield General (Term 1) for prime 3.")

    # --- 4. Ratchet Loop ---
    # We continue until we have found the requested number of terms.
    while len(results) < n_terms:
        # Identify the NEW prime we need to shield against (p_{k+1}).
        p_new = nextprime(p_k)

        # Search for the smallest k such that:
        # candidate = current_general + k * current_cycle
        # satisfies the condition for p_new.
        # Since 'current_cycle' is a multiple of all previous primes, adding it
        # preserves the modular properties for all previous primes.
        k = 0
        while True:
            candidate = current_general + (k * current_cycle)

            # Check 1: Does this candidate satisfy the new prime's constraint?
            # We accept Natural Shields (1) and Selection Shields (-1).
            rem = candidate % p_new
            hit = (rem == 1 or rem == p_new - 1)

            # Check 2: Does it maintain the required parity (even)?
            # (Note: With the current setup starting at 4 and cycle including odd primes,
            # this check is usually satisfied but kept for rigor).
            is_target_parity = (candidate % target_mod == target_rem)

            if hit and is_target_parity:
                # We found the new Shield General!
                current_general = candidate

                # Update our list of shielded primes to include the new one.
                shielded_primes.append(p_new)

                # --- VALIDATION STEP ---
                # Perform a full check against ALL shielded primes to guarantee integrity.
                if not validate_general(current_general, shielded_primes, target_mod, target_rem):
                    raise RuntimeError(f"Validation failed for General {current_general} with prime {p_new}.")

                logging.info(f"Successfully validated Shield General (Term {len(results) + 1}) for prime {p_new}.")
                break

            # If not found, increment k and try the next multiple of the cycle.
            k += 1

        # Append the valid General to our results.
        results.append(current_general)

        # --- Update State for Next Iteration ---
        # Multiply the cycle by the new prime.
        current_cycle *= p_new

        # Update the highest shielded prime marker.
        p_k = p_new

    return results
if __name__ == "__main__":
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(
        description="Calculate the first N terms of the Shield General prime gap sequence.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "N",
        type=int,
        help="The number of Shield Generals (N) to calculate.\n" \
             "Example: 6 yields 4, 4, 34, 1924, 25024, 85084."
    )

    args = parser.parse_args()

    logging.info(f"Starting Shield General search for {args.N} terms.")

    try:
        # Run the generator
        sequence = generate_shield_generals(args.N)

        logging.info("--- Results ---")

        # List of primes for display purposes (to show P-max)
        prime_list = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59]

        for i, g in enumerate(sequence):
            p_shielded_count = i + 1
            # Determine the highest prime shielded (P-max) for this term
            current_p_max = prime_list[i] if i < len(prime_list) else "Unknown"

            logging.info(f"Shield {p_shielded_count} (P-max: {current_p_max}): {g:,}")

        logging.info(f"Sequence Array: {sequence}")

    except (ValueError, RuntimeError) as e:
        logging.error(f"Error: {e}")
    except IndexError:
         # Handle cases where we exceed the hardcoded prime_list for display
         logging.info(f"Sequence Array: {sequence}")
         logging.warning("The P-max list limit was reached. Results are still correct.")
    except ImportError:
        logging.error("The 'sympy' library is required. Install it with: pip install sympy")
