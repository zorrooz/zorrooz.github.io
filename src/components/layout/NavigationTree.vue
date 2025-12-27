<!-- NavigationTree.vue -->
<template>
  <div class="navigation-tree">
    <div v-for="category in navigationTree" :key="category.name" class="category-group">
      <h3 class="category-name">{{ category.name }}</h3>
      <ul class="article-list article-list-root">
        <li v-if="category.children && category.children.length" v-for="dir in category.children" :key="dir.name"
          class="article-item">
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
        <li v-if="category.files && category.files.length" v-for="file in category.files" :key="file.path"
          class="article-item">
          <router-link :to="toArticle(file.path)" class="article-link level-3" :class="{ active: isActive(file.path) }">
            {{ file.title }}
          </router-link>
        </li>
      </ul>
    </div>
  </div>
</template>

<script>
import { useI18n } from 'vue-i18n'
import { loadCategories } from '@/utils/contentLoader';

/* 
  NavigationTree
  - 根据 categories 数据构建分级导航树
*/
export default {
  name: 'NavigationTree',
  setup() {
    const { locale } = useI18n()
    return { locale }
  },
  data() {
    return {
      navigationTree: [],
      currentPath: '',
      categoryData: []
    };
  },
  watch: {
    '$route': {
      handler(to) {
        this.currentPath = to.params.path ? to.params.path.join('/') : '';
        this.buildTree();
      },
      immediate: true,
      deep: true
    },
    locale() { this.loadCategoryData().then(() => { this.currentPath = this.$route.params.path ? this.$route.params.path.join('/') : ''; this.buildTree(); }); }
  },
  async mounted() {
    await this.loadCategoryData();
    this.buildTree();
  },
  methods: {
    async loadCategoryData() {
      try {
        this.categoryData = await loadCategories() || [];
      } catch (error) {
        console.error('Failed to load category data:', error);
        this.categoryData = [];
      }
      return Promise.resolve();
    },
    toArticle(path) {
      return { name: 'Article', params: { path: path.replace(/\.md$/, '').split('/') } };
    },
    isActive(path) {
      const currentPath = this.currentPath.replace(/\.md$/, '');
      const articlePath = path.replace(/\.md$/, '');

      const cleanCurrentPath = currentPath.replace(/-en$/, '');
      const cleanArticlePath = articlePath.replace(/-en$/, '');

      return cleanCurrentPath === cleanArticlePath;
    },
    buildTree() {
      const path = this.currentPath; if (!path) { this.navigationTree = []; return; }
      const segs = path.split('/').filter(Boolean); if (segs.length < 2) { this.navigationTree = []; return; }
      const type = segs[0]; const group = segs[1];

      let targetItem = null;
      if (Array.isArray(this.categoryData)) {
        outer:
        for (const section of this.categoryData) {
          if (!Array.isArray(section.items)) continue;
          for (const item of section.items) {
            if (item?.name !== group) continue;

            let hasTypeMatch = false;

            if (Array.isArray(item.articles)) {
              hasTypeMatch = item.articles.some(a => typeof a?.articleUrl === 'string' && a.articleUrl.includes(`/article/${type}/`));
            }
            if (!hasTypeMatch && Array.isArray(item.categories)) {
              hasTypeMatch = item.categories.some(cat =>
                Array.isArray(cat.articles) &&
                cat.articles.some(a => typeof a?.articleUrl === 'string' && a.articleUrl.includes(`/article/${type}/${group}/`))
              );
            }

            if (hasTypeMatch) {
              targetItem = item;
              break outer;
            }
          }
        }
      }
      if (!targetItem) { this.navigationTree = []; return; }

      const rootFiles = [];
      const children = [];

      const toFile = (title, articleUrl) => {
        const parts = String(articleUrl).replace(/^\/+/, '').split('/');
        const i0 = parts[0] === 'article' ? 1 : 0;
        const t = parts[i0];
        const g = parts[i0 + 1];
        const rest = parts.slice(i0 + 2);
        const pathNoExt = `${t}/${g}/${rest.join('/')}`;
        return { title, path: `${pathNoExt}.md` };
      };

      if (Array.isArray(targetItem.articles)) {
        targetItem.articles.forEach(a => { if (a?.articleUrl) rootFiles.push(toFile(a.title, a.articleUrl)); });
      }

      if (Array.isArray(targetItem.categories)) {
        targetItem.categories.forEach(cat => {
          const files = [];
          if (Array.isArray(cat.articles)) {
            cat.articles.forEach(a => { if (a?.articleUrl) files.push(toFile(a.title, a.articleUrl)); });
          }
          if (files.length) {
            children.push({
              name: cat.title || cat.key,
              type: 'directory',
              files
            });
          }
        });
      }

      this.navigationTree = [{
        name: targetItem.title || targetItem.name || group,
        files: rootFiles,
        children
      }];
    }
  }
};
</script>

<style scoped>
.navigation-tree {
  padding: 0;
  font-size: 0.95rem;
}

.category-group {
  margin-bottom: 1rem;
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
