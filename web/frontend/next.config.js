/** @type {import('next').NextConfig} */
const BACKEND = process.env.VAXFORGE_API || "http://127.0.0.1:8011";

const nextConfig = {
  reactStrictMode: true,
  // Frontend'ten /api/* çağrıları backend'e proxy'lenir (aynı origin, CORS derdi yok).
  async rewrites() {
    return [{ source: "/api/:path*", destination: `${BACKEND}/api/:path*` }];
  },
};

module.exports = nextConfig;
