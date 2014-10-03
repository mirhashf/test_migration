JAX-RS API mostly intended for interviews etc.
---------------------------------------------
1. /api/mutations - Returns mutations. Can be paged by appending ?page=1&pagesize=20

2. /api/mutations/{id} - Returns a single mutation by id

3. /api/mutations/{id}/annotations - Returns annotations for a specific mutation

4. /api/test - Health check

To run the API, build with Maven and run the produced WAR file in the target directory, or run **com.bina.mutations.app.Runner.java** from an IDE like Eclipse.

**Note:** This API is open to for cross-origin resource sharing (CORS) so AJAX calls from an HTML5 compliant browser can call this API.

Possible interview questions
----------------------------
1. Using the given api, develop a UI that prints the first 10 mutations and their respective annotations 

1a. Print the mutations in chromosome-position ascending order

1b. Filter out all annotations that are TOLERATED