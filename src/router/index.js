// @ts-check
import HomeView from '@/views/Home.vue'
import CategoryView from '@/views/Category.vue'
import ResourceView from '@/views/Resource.vue'
import AboutView from '@/views/About.vue'
import ArticleView from '@/views/Article.vue'

/** Resolve the preferred locale prefix (client-side only; SSR defaults to /zh). */
const localePrefix = () => {
  if (typeof window !== 'undefined') {
    const saved = localStorage.getItem('locale')
    if (saved === 'en-US') return '/en'
    if (saved === 'zh-CN') return '/zh'
    if (navigator.language && navigator.language.toLowerCase().startsWith('en')) return '/en'
  }
  return '/zh'
}

/** @param {string} prefix */
const prefixedRoutes = (prefix) => [
  { path: `${prefix}/`, name: `${prefix}-Home`, component: HomeView },
  { path: `${prefix}/category`, name: `${prefix}-Category`, component: CategoryView },
  { path: `${prefix}/resource`, name: `${prefix}-Resource`, component: ResourceView },
  { path: `${prefix}/about`, name: `${prefix}-About`, component: AboutView },
  { path: `${prefix}/article/:path*`, name: `${prefix}-Article`, component: ArticleView, props: true },
]

/** @type {import('vue-router').RouteRecordRaw[]} */
export const routes = [
  // Legacy unprefixed URLs redirect to the locale-prefixed equivalent
  { path: '/', redirect: () => localePrefix() },
  { path: '/category', redirect: (to) => `${localePrefix()}${to.path}` },
  { path: '/resource', redirect: (to) => `${localePrefix()}${to.path}` },
  { path: '/about', redirect: (to) => `${localePrefix()}${to.path}` },
  { path: '/article/:path*', redirect: (to) => `${localePrefix()}${to.path}` },
  ...prefixedRoutes('/zh'),
  ...prefixedRoutes('/en'),
]
