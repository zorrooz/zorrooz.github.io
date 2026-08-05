import HomeView from '@/views/Home.vue'
import CategoryView from '@/views/Category.vue'
import ResourceView from '@/views/Resource.vue'
import AboutView from '@/views/About.vue'
import ArticleView from '@/views/Article.vue'

import type { RouteRecordRaw } from 'vue-router'

import { LOCALE_SEGMENTS, preferredLocaleSegment, type LocaleSegment } from '@/config/site'

const prefixedRoutes = (prefix: LocaleSegment): RouteRecordRaw[] => [
  { path: `/${prefix}/`, name: `${prefix}-Home`, component: HomeView },
  { path: `/${prefix}/category`, name: `${prefix}-Category`, component: CategoryView },
  { path: `/${prefix}/resource`, name: `${prefix}-Resource`, component: ResourceView },
  { path: `/${prefix}/about`, name: `${prefix}-About`, component: AboutView },
  {
    path: `/${prefix}/article/:path*`,
    name: `${prefix}-Article`,
    component: ArticleView,
    props: true,
  },
]

export const routes: RouteRecordRaw[] = [
  // Legacy unprefixed URLs redirect to the locale-prefixed equivalent
  { path: '/', redirect: () => `/${preferredLocaleSegment()}/` },
  { path: '/category', redirect: (to) => `/${preferredLocaleSegment()}${to.path}` },
  { path: '/resource', redirect: (to) => `/${preferredLocaleSegment()}${to.path}` },
  { path: '/about', redirect: (to) => `/${preferredLocaleSegment()}${to.path}` },
  { path: '/article/:path*', redirect: (to) => `/${preferredLocaleSegment()}${to.path}` },
  ...LOCALE_SEGMENTS.map(prefixedRoutes).flat(),
]
