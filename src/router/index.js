import { createRouter, createWebHashHistory } from 'vue-router'
import HomeView from '@/views/Home.vue'
import CategoryView from '@/views/Category.vue'
import ResourceView from '@/views/Resource.vue'
import AboutView from '@/views/About.vue'
import ArticleView from '@/views/Article.vue'

const router = createRouter({
  history: createWebHashHistory(import.meta.env.BASE_URL),
  routes: [{
    path: '/',
    name: 'Home',
    component: HomeView
  }, {
    path: '/category',
    name: 'Category',
    component: CategoryView
  }, {
    path: '/resource',
    name: 'Resource',
    component: ResourceView
  }, {
    path: '/about',
    name: 'About',
    component: AboutView
  }, {
    path: '/article',
    name: 'Article',
    component: ArticleView
  }
  ],
})

export default router
