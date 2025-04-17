# --- Finding the Angle at which the Perimeter of the Ellipse is Divided Equally ---
# This function computes the angle theta where the arc length from (-a, 0) to a point on the ellipse
# equals 1/4 of the total perimeter, using numerical integration and the bisection method.
# Inputs: a (semi-major axis), b (semi-minor axis)
# Output: theta (angle in radians)
function perimeter_angle_find(a::Float64, b::Float64)
    # ---- Integral Function ----
    # Define a helper function to compute the arc length of the ellipse from -lim to lim.
    # The arc length integral is derived from the ellipse equation x^2/a^2 + y^2/b^2 = 1.
    # Parameters:
    #   k: semi-major axis (a)
    #   l: semi-minor axis (b)
    #   lim: upper limit of integration (x-coordinate from -lim to lim)
    function integral_curve_length(k, l, lim)
        # Use quadgk for numerical integration of the arc length formula:
        # ds = sqrt(1 + (dy/dx)^2) dx, where dy/dx is derived from the ellipse equation.
        # For ellipse x^2/k^2 + y^2/l^2 = 1, dy/dx = -(l^2 x) / (k^2 y), and y = l * sqrt(1 - x^2/k^2).
        # The integrand becomes sqrt(1 + x^2 * l^2 / (k^4 - k^2 * x^2)).
        out, _ = quadgk(x -> sqrt(1 + x^2 * l^2 / (k^4 - k^2 * x^2)), -lim, lim, rtol=1e-8)
        return out  # Return the computed arc length
    end

    # ---- Compute Full Perimeter ----
    # The total perimeter of the ellipse is twice the arc length from -a to a (due to symmetry).
    peri = 2 * integral_curve_length(a, b, a)
    # Note: This approximates the full perimeter using numerical integration.

    # ---- Define Objective Function ----
    # Define a function to find the difference between the arc length from -lim to lim
    # and the target arc length (1/4 of the perimeter).
    # Parameters:
    #   lim: x-coordinate limit to test
    #   k, l: semi-major and semi-minor axes (a, b)
    #   final: target arc length (peri/4)
    function f(lim, k, l, final)
        val = integral_curve_length(k, l, lim) - final
        return val  # Positive if arc length > target, negative if < target
    end

    # ---- Bisection Method ----
    # Goal: Find the x-coordinate (a_mid) where the arc length from -a_mid to a_mid
    # equals 1/4 of the total perimeter, then compute the corresponding angle.
    target = peri / 4.0  # Target arc length is 1/4 of the total perimeter

    # ---- Initial Bounds ----
    # Set initial bounds for the bisection method:
    a_lower = 0.0    # Lower bound (x = 0, start of ellipse at (0, b))
    a_upper = a      # Upper bound (x = a, end of first quadrant at (a, 0))

    # ---- Verify Root Existence ----
    # Evaluate the function at the bounds to ensure a root exists (sign change).
    f_lower = f(a_lower, a, b, target)  # Arc length at x = 0 (should be 0 - target < 0)
    f_upper = f(a_upper, a, b, target)  # Arc length at x = a (should be peri/2 - target > 0)

    # Check if the function changes sign between bounds, a requirement for bisection.
    if f_lower * f_upper >= 0
        error("Function values at bounds have the same sign. Adjust bounds.")
        # This ensures the target lies between f_lower and f_upper.
    end

    # ---- Bisection Parameters ----
    tolerance = 1e-6       # Convergence tolerance (stop when solution is within 1e-6)
    max_iterations = 1000  # Maximum iterations to prevent infinite loops
    a_mid = (a_lower + a_upper) / 2.0  # Initial guess: midpoint of bounds

    # ---- Bisection Loop ----
    # Iteratively narrow the bounds to find the x-coordinate where arc length = target.
    for i in 1:max_iterations
        f_mid = f(a_mid, a, b, target)  # Evaluate function at current midpoint
        
        # Check convergence: stop if f_mid is near zero or bounds are sufficiently close
        if abs(f_mid) < tolerance || (a_upper - a_lower) / 2.0 < tolerance
            break  # Exit loop when solution is found
        end
        
        # Update bounds based on the sign of f_mid:
        if f_mid > 0
            a_upper = a_mid  # Arc length too large, move upper bound down
        else
            a_lower = a_mid  # Arc length too small, move lower bound up
        end
        a_mid = (a_lower + a_upper) / 2.0  # Compute new midpoint
    end

    # ---- Compute Corresponding y and Angle ----
    # Given x = a_mid, compute y using the ellipse equation: y = b * sqrt(1 - x^2/a^2)
    b_mid = b * sqrt(1 - a_mid^2 / a^2)  # y-coordinate at x = a_mid (positive quadrant)
    
    # Compute the angle theta using arctangent: theta = atan(y/x)
    theta = atan(b_mid / a_mid)  # Angle in radians from (0,0) to (a_mid, b_mid)

    return theta  # Return the angle where the perimeter is divided into quarters
end