<!-- NavigationTree.vue -->
<template>
  <div class="navigation-tree">
    <div v-for="category in navigationTree" :key="category.name" class="tree-group">
      <div class="tree-label">{{ category.name }}</div>
      <ul class="article-list article-list-root">
        <template v-for="dir in category.children" :key="dir.name">
          <li v-if="dir.files && dir.files.length" class="article-item">
            <div class="directory-node">
              <div class="tree-item tree-item--folder">
                {{ dir.name }}
              </div>
              <ul v-if="dir.files && dir.files.length" class="article-list sub-list files-level">
                <li v-for="file in dir.files" :key="file.path" class="article-item">
                  <router-link :to="toArticle(file.path)" class="tree-item tree-item--child"
                    :class="{ 'tree-item--active': isActive(file.path) }">
                    {{ file.title }}
                  </router-link>
                </li>
              </ul>
              <ul v-if="dir.children && dir.children.length" class="article-list sub-list">
                <li v-for="sub in dir.children" :key="sub.name" class="article-item">
                  <div class="directory-node">
                    <div class="tree-item tree-item--folder">
                      {{ sub.name }}
                    </div>
                    <ul v-if="sub.files && sub.files.length" class="article-list sub-list files-level">
                      <li v-for="file in sub.files" :key="file.path" class="article-item">
                        <router-link :to="toArticle(file.path)" class="tree-item tree-item--child"
                          :class="{ 'tree-item--active': isActive(file.path) }">
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
          <router-link :to="toArticle(file.path)" class="tree-item tree-item--child"
            :class="{ 'tree-item--active': isActive(file.path) }">
            {{ file.title }}
          </router-link>
        </li>
      </ul>
    </div>
  </div>
</template>

<script setup lang="ts">
import { onMounted, ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { useRoute } from 'vue-router'
import { loadCategories } from '@/utils/contentLoader'
import { toLocalePath } from '@/utils/localePath'

interface TreeFile {
  title: string
  path: string
}

interface TreeDir {
  name: string
  type: string
  files: TreeFile[]
  children?: TreeDir[]
}

interface TreeCategory {
  name: string
  files: TreeFile[]
  children: TreeDir[]
}

const { locale } = useI18n()
const route = useRoute()

const navigationTree = ref<TreeCategory[]>([])
const currentPath = ref('')
const categoryData = ref<any[]>([])

function loadCategoryData() {
  try {
    categoryData.value = loadCategories() || []
  } catch (error) {
    console.error('Failed to load category data:', error)
    categoryData.value = []
  }
}

function toArticle(path: string) {
  return toLocalePath(`/article/${path.replace(/\.md$/, '')}`)
}

function isActive(path: string) {
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

  let targetItem: any = null
  if (Array.isArray(categoryData.value)) {
    outer:
    for (const section of categoryData.value) {
      if (!Array.isArray(section.items)) continue
      for (const item of section.items) {
        if (item?.name !== group) continue

        let hasTypeMatch = false

        if (Array.isArray(item.articles)) {
          hasTypeMatch = item.articles.some((a: any) => typeof a?.articleUrl === 'string' && a.articleUrl.includes(`/article/${type}/`))
        }
        if (!hasTypeMatch && Array.isArray(item.categories)) {
          hasTypeMatch = item.categories.some((cat: any) =>
            Array.isArray(cat.articles) &&
            cat.articles.some((a: any) => typeof a?.articleUrl === 'string' && a.articleUrl.includes(`/article/${type}/${group}/`))
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

  const rootFiles: TreeFile[] = []
  const children: TreeDir[] = []

  const toFile = (title: string, articleUrl: string): TreeFile => {
    const parts = String(articleUrl).replace(/^\/+/, '').split('/')
    const i0 = parts[0] === 'article' ? 1 : 0
    const t = parts[i0]
    const g = parts[i0 + 1]
    const rest = parts.slice(i0 + 2)
    const pathNoExt = `${t}/${g}/${rest.join('/')}`
    return { title, path: `${pathNoExt}.md` }
  }

  if (Array.isArray(targetItem.articles)) {
    targetItem.articles.forEach((a: any) => { if (a?.articleUrl) rootFiles.push(toFile(a.title, a.articleUrl)) })
  }

  if (Array.isArray(targetItem.categories)) {
    targetItem.categories.forEach((cat: any) => {
      const files: TreeFile[] = []
      if (Array.isArray(cat.articles)) {
        cat.articles.forEach((a: any) => { if (a?.articleUrl) files.push(toFile(a.title, a.articleUrl)) })
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
  const p = route.params.path
  currentPath.value = Array.isArray(p) ? p.join('/') : (p || '')
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
  padding: 0.25rem 0;
  font-size: var(--text-base);
}

.tree-group {
  margin-bottom: var(--sp-6);
}

.tree-label {
  font-size: 12px;
  font-weight: 600;
  letter-spacing: 0.08em;
  text-transform: uppercase;
  color: var(--fg-3);
  padding: 0 var(--sp-2);
  margin-bottom: var(--sp-2);
  display: flex;
  align-items: center;
  gap: var(--sp-2);
}

.tree-label::after {
  content: '';
  flex: 1;
  height: 1px;
  background: var(--line);
}

.article-list {
  list-style: none;
  padding-left: 0;
  margin: 0;
}

.article-list-root {
  padding-left: 0;
}

.article-item {
  margin-bottom: 1px;
}

.directory-node {
  padding: 0.1rem 0;
}

.tree-item {
  display: flex;
  align-items: center;
  gap: 8px;
  font-size: var(--text-sm);
  color: var(--fg-2);
  padding: 7px 10px;
  border-radius: var(--radius-sm);
  margin-bottom: 1px;
  transition: color 0.14s ease, background-color 0.14s ease;
  cursor: pointer;
  width: 100%;
  text-align: left;
  border: none;
  background: transparent;
  position: relative;
  text-decoration: none;
  line-height: 1.45;
}

.tree-item:hover {
  background: var(--surface-2);
  color: var(--fg);
}

/* 当前文章 hover 保持蓝色，不落入 hover 的黑色 */
.tree-item--active:hover {
  color: var(--primary);
  background: var(--tint-strong);
}

.tree-item--active {
  color: var(--primary);
  font-weight: 600;
  background: var(--tint);
}

.tree-item--folder {
  font-weight: 600;
  color: var(--fg);
  cursor: default;
}

.tree-item--child {
  padding-left: 10px;
  font-size: var(--text-sm);
}

/* 子目录内的文章（二级及更深）比分类直属文章小一档 */
.files-level .tree-item--child {
  font-size: 13px;
}

.sub-list {
  padding-left: 0;
}
</style>
