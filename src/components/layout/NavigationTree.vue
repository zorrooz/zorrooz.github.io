<!-- NavigationTree.vue -->
<template>
  <div class="navigation-tree">
    <div v-for="category in navigationTree" :key="category.name" class="category-group">
      <h3 class="category-name">{{ category.name }}</h3>
      <ul class="article-list article-list-root">
        <template v-for="dir in category.children" :key="dir.name">
          <li v-if="dir.files && dir.files.length" class="article-item">
            <div class="directory-node">
              <span class="directory-name level-2">{{ dir.name }}</span>
              <ul v-if="dir.files && dir.files.length" class="article-list sub-list files-level">
                <li v-for="file in dir.files" :key="file.path" class="article-item">
                  <router-link :to="toArticle(file.path)" class="article-link level-3"
                    :class="{ active: isActive(file.path) }">
                    {{ file.title }}
                  </router-link>
                </li>
              </ul>
              <ul v-if="dir.children && dir.children.length" class="article-list sub-list">
                <li v-for="sub in dir.children" :key="sub.name" class="article-item">
                  <div class="directory-node">
                    <span class="directory-name level-2">{{ sub.name }}</span>
                    <ul v-if="sub.files && sub.files.length" class="article-list sub-list files-level">
                      <li v-for="file in sub.files" :key="file.path" class="article-item">
                        <router-link :to="toArticle(file.path)" class="article-link level-3"
                          :class="{ active: isActive(file.path) }">
                          {{ file.title }}
                        </router-link>
                      </li>
                    </ul>
                  </div>
                </li>
              </ul>
            </div>
          </li>
        </template>
        <li v-for="file in category.files" :key="file.path" class="article-item">
          <router-link :to="toArticle(file.path)" class="article-link level-3" :class="{ active: isActive(file.path) }">
            {{ file.title }}
          </router-link>
        </li>
      </ul>
    </div>
  </div>
</template>

<script setup>
import { onMounted, ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { useRoute } from 'vue-router'
import { loadCategories } from '@/utils/contentLoader'

const { locale } = useI18n()
const route = useRoute()

const navigationTree = ref([])
const currentPath = ref('')
const categoryData = ref([])

function loadCategoryData() {
  try {
    categoryData.value = loadCategories() || []
  } catch (error) {
    console.error('Failed to load category data:', error)
    categoryData.value = []
  }
}

function toArticle(path) {
  return { name: 'Article', params: { path: path.replace(/\.md$/, '').split('/') } }
}

function isActive(path) {
  const current = currentPath.value.replace(/\.md$/, '')
  const articlePath = path.replace(/\.md$/, '')

  const cleanCurrentPath = current.replace(/-en$/, '')
  const cleanArticlePath = articlePath.replace(/-en$/, '')

  return cleanCurrentPath === cleanArticlePath
}

function buildTree() {
  const path = currentPath.value
  if (!path) { navigationTree.value = []; return }
  const segs = path.split('/').filter(Boolean)
  if (segs.length < 2) { navigationTree.value = []; return }
  const type = segs[0]
  const group = segs[1]

  let targetItem = null
  if (Array.isArray(categoryData.value)) {
    outer:
    for (const section of categoryData.value) {
      if (!Array.isArray(section.items)) continue
      for (const item of section.items) {
        if (item?.name !== group) continue

        let hasTypeMatch = false

        if (Array.isArray(item.articles)) {
          hasTypeMatch = item.articles.some(a => typeof a?.articleUrl === 'string' && a.articleUrl.includes(`/article/${type}/`))
        }
        if (!hasTypeMatch && Array.isArray(item.categories)) {
          hasTypeMatch = item.categories.some(cat =>
            Array.isArray(cat.articles) &&
            cat.articles.some(a => typeof a?.articleUrl === 'string' && a.articleUrl.includes(`/article/${type}/${group}/`))
          )
        }

        if (hasTypeMatch) {
          targetItem = item
          break outer
        }
      }
    }
  }
  if (!targetItem) { navigationTree.value = []; return }

  const rootFiles = []
  const children = []

  const toFile = (title, articleUrl) => {
    const parts = String(articleUrl).replace(/^\/+/, '').split('/')
    const i0 = parts[0] === 'article' ? 1 : 0
    const t = parts[i0]
    const g = parts[i0 + 1]
    const rest = parts.slice(i0 + 2)
    const pathNoExt = `${t}/${g}/${rest.join('/')}`
    return { title, path: `${pathNoExt}.md` }
  }

  if (Array.isArray(targetItem.articles)) {
    targetItem.articles.forEach(a => { if (a?.articleUrl) rootFiles.push(toFile(a.title, a.articleUrl)) })
  }

  if (Array.isArray(targetItem.categories)) {
    targetItem.categories.forEach(cat => {
      const files = []
      if (Array.isArray(cat.articles)) {
        cat.articles.forEach(a => { if (a?.articleUrl) files.push(toFile(a.title, a.articleUrl)) })
      }
      if (files.length) {
        children.push({
          name: cat.title || cat.key,
          type: 'directory',
          files
        })
      }
    })
  }

  navigationTree.value = [{
    name: targetItem.title || targetItem.name || group,
    files: rootFiles,
    children
  }]
}

function syncPathAndBuild() {
  currentPath.value = route.params.path ? route.params.path.join('/') : ''
  buildTree()
}

watch(() => route.params.path, syncPathAndBuild, { immediate: true })

watch(locale, () => {
  loadCategoryData()
  syncPathAndBuild()
})

onMounted(() => {
  loadCategoryData()
  buildTree()
})
</script>

<style scoped>
.navigation-tree {
  padding: 0.5rem 0;
  font-size: 0.95rem;
}

.category-group {
  margin-bottom: 0.5rem;
}

.category-group:last-child {
  margin-bottom: 0;
}

.category-name {
  font-size: 1.1rem;
  font-weight: 700;
  color: var(--app-text-muted);
  margin-bottom: 0.5rem;
  padding-bottom: 0.25rem;
  padding-left: 0.75rem;
}

.article-list {
  list-style: none;
  padding-left: 0;
}

.article-list-root {
  padding-left: 0;
}

.directory-node {
  padding: 0.25rem 0;
}

.directory-name.level-2 {
  display: block;
  font-weight: 700;
  color: var(--app-text-emphasis);
  margin-left: 0;
  padding-left: 0.75rem;
}

.sub-list {
  padding-left: 0.75rem;
  margin-top: 0.5rem;
}

.files-level {
  padding-left: 0.75rem;
}

.article-item {
  margin-bottom: 0.3rem;
}

.article-link {
  display: block;
  padding: 0.5rem 0.75rem;
  text-decoration: none;
  color: var(--app-text-muted);
  border-radius: 0.25rem;
  transition: background-color 0.2s ease, color 0.2s ease;
}

.article-link:hover,
.article-link:focus {
  background-color: var(--app-primary-bg-subtle);
  color: var(--app-primary);
}

.article-link.active {
  background-color: transparent;
  color: var(--app-primary);
  font-weight: 700;
}
</style>
