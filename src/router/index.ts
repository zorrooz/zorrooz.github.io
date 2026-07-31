import HomeView from '@/views/Home.vue'
import CategoryView from '@/views/Category.vue'
import ResourceView from '@/views/Resource.vue'
import AboutView from '@/views/About.vue'
import ArticleView from '@/views/Article.vue'

import type { RouteRecordRaw } from 'vue-router'

import { LOCALE_PREFIXES, SUPPORTED_LOCALES } from '@/config/site.ts'

/** Resolve the preferred locale prefix (client-side only; SSR defaults to /zh). */
const localePrefix = (): string => {
  if (typeof window !== 'undefined') {
    const saved = localStorage.getItem('locale')
    if (saved === SUPPORTED_LOCALES[1]) return LOCALE_PREFIXES[1]
    if (saved === SUPPORTED_LOCALES[0]) return LOCALE_PREFIXES[0]
    if (navigator.language && navigator.language.toLowerCase().startsWith('en')) return LOCALE_PREFIXES[1]
  }
  return LOCALE_PREFIXES[0]
}

const prefixedRoutes = (prefix: (typeof LOCALE_PREFIXES)[number]): RouteRecordRaw[] => [
  { path: `${prefix}/`, name: `${prefix}-Home`, component: HomeView },
  { path: `${prefix}/category`, name: `${prefix}-Category`, component: CategoryView },
  { path: `${prefix}/resource`, name: `${prefix}-Resource`, component: ResourceView },
  { path: `${prefix}/about`, name: `${prefix}-About`, component: AboutView },
  { path: `${prefix}/article/:path*`, name: `${prefix}-Article`, component: ArticleView, props: true },
]

export const routes: RouteRecordRaw[] = [
  // Legacy unprefixed URLs redirect to the locale-prefixed equivalent
  { path: '/', redirect: () => localePrefix() },
  { path: '/category', redirect: (to) => `${localePrefix()}${to.path}` },
  { path: '/resource', redirect: (to) => `${localePrefix()}${to.path}` },
  { path: '/about', redirect: (to) => `${localePrefix()}${to.path}` },
  { path: '/article/:path*', redirect: (to) => `${localePrefix()}${to.path}` },
  ...prefixedRoutes(LOCALE_PREFIXES[0]),
  ...prefixedRoutes(LOCALE_PREFIXES[1]),
]
