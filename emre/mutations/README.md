JAX-RS API mostly intended for interviews etc.
---------------------------------------------
1. /api/mutations - Returns mutations. Can be paged by appending ?page=1&pagesize=20

2. /api/mutations/{id}/annotations - Returns annotations for a specific mutation

3. /api/test - Health check

**Note:** This API is open to for cross-origin resource sharing (CORS) so AJAX calls from an HTML5 compliant browser can call this API.